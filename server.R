library(latticeExtra)
library(Kmisc)
library(rCharts)
library(flowWorkspace)
library(gridExtra)
library(hexbin)
library(ggplot2)
library(data.table.extras)

## debug
if (FALSE) {
  channel_1 <- "Live"
  channel_2 <- "Dead"
  GS <- gs <- load_gs("gs/LoveLab-Gated")
  FS <- fs <- getData(gs)
}

options(debug=TRUE)
debug <- function(fmt, ...) {
  if (isTRUE(getOption("debug"))) {
    sprintf(fmt, ...)
  }
  return( invisible(NULL) )
}

file.ext <- function(x) {
  return( gsub(".*\\.", "", x) )
}

shinyServer( function(input, output, session) {
  
  GS <- gs <- NULL
  FS <- fs <- NULL
  template <- NULL
  
  ## Utility Functions
  getSampleSelection <- reactive({
    samples <- input$sample_selection
    
    if (!length(samples))
      samples <- sampleNames(gs)
    
    return(samples)
  })
  
  ## An observer for the 'select all samples' button.
  observe({
    if (input$select_all_samples == 0) {
      return( invisible(NULL) )
    }
    
    isolate({
      updateSelectInput(session, "sample_selection", choices=sort(sampleNames(gs)), selected='')
    })
    
  })
  
  ## An observer for the remove gate button.
  observe({
    debug("Removing gate!")
    if (input$remove_gate == 0) 
      return( invisible(NULL) )
    
    isolate({
      
      ## remove the node from the gs
      parent <- input$parent_selection
      samples <- getSampleSelection()
      
      ## only remove gates for samples that have gates at this node
      missing_samples <- NULL
      
      for (sample in samples) {
        if (!(parent %in% getNodes(gs[[sample]], isPath=TRUE))) {
          missing_samples <- c(missing_samples, sample)
        }
      }
      
      samples <- samples[ !(samples %in% missing_samples) ]
      
      pp <- getParent(gs[[samples[1]]], parent, isPath=TRUE)
      Rm(parent, gs[samples])
      
      ## remove the node from the template
      
      ## algorithm
      ## check the template and get all entries that have 'sample'
      ## as an entry in the samples column. if the parent is equal to 'parent',
      ## remove sample from that entry.
      ## at the end, if the sample column is empty, delete it
      template_samples <- strsplit( as.character(template$samples), ",")
      for (sample in samples) {
        for (i in 1:nrow(template)) {
          if (template$parent[i] == pp && sample %in% template_samples[[i]]) {
            template_samples[[i]] <- template_samples[[i]][ template_samples[[i]] != sample ]
          }
        }
      }
      template$samples <- sapply(template_samples, function(x) {
        paste(x, collapse=",")
      })
      template <<- template[ template$samples != '', ]
      
      ## update the parent selection
      allNodes <- unique( unlist( lapply(gs, function(x) getNodes(x, isPath=TRUE) ) ) )
      updateSelectInput(session, "parent_selection", choices=allNodes, selected=pp)
      
    })
    
  })
  
  observe({
    debug("Trying to read and update select inputs.")
    file <- input$data_selection
    GS <<- gs <<- load_gs(file)
    FS <<- fs <<- getData(gs)
    template <<- vector("list", 11)
    names(template) <<- c("alias", "pop", "parent", "dims", "gating_method", 
      "gating_args", "collapseDataForGating", "groupBy", "preprocessing_method", 
      "preprocessing_args", "samples")
    template <<- as.data.frame(template)
    channels <- sort(getData(gs)@colnames)
    nodes <- getNodes(gs[[1]], isPath=TRUE)
    updateSelectInput(session, "parent_selection", choices=nodes)
    updateSelectInput(session, "channel_1", choices=channels)
    updateSelectInput(session, "channel_2", choices=channels)
    updateSelectInput(session, "sample_selection", choices=sort(sampleNames(gs)))
  })
  
  ## An observer for the button that applies a gate.
  observe({
    debug("Trying to apply gate.")
    if (input$apply_gate == 0) 
      return(NULL)
    
    isolate({
      channel_1 <- input$channel_1
      channel_2 <- input$channel_2
      parent <- input$parent_selection
      gating_method <- input$gating_method
      gate_name <- input$gate_name
      samples <- getSampleSelection()
      selected_channel <- input$select_channels
      gate_both <- input$gate_positive_and_negative
      
      tryCatch({
        fn <- get(gating_method, envir=asNamespace("openCyto"))
        fnArgs <- formals(fn)
        if (gating_method == "flowClust.2d") {
          fnArgs <- fnArgs[5:length(fnArgs)]
        } else {
          fnArgs <- fnArgs[4:length(fnArgs)]
        }
        fnArgs <- fnArgs[ names(fnArgs) != "..." ]
        args <- vector("list", length(fnArgs))
        names(args) <- names(fnArgs)
        for (i in seq_along(fnArgs)) {
          tmp <- eval(parse(text=input[[paste(gating_method, names(fnArgs)[i], sep="_")]]))
          if (is.null(tmp)) {
            args[i] <- list(NULL)
          } else {
            args[[i]] <- tmp
          }
        }
        
        if (gating_method == "flowClust.1d" && is.null(args[["K"]])) {
          warning("Changing 'K' from NULL to 2")
          args[["K"]] <- 2
        }
        
        ## the first three arguments are fixed and represent the flow frame,
        ## the channel name, and the filterId, for the 1d gates
        
        if (gating_method == "flowClust.2d") {
          args <- c('', '', '', '', args)
          args[[2]] <- channel_1
          args[[3]] <- channel_2
          args[[4]] <- paste(gating_method, paste(channel_1, channel_2, sep=":"), sep="_")
        } else {
          args <- c('', '', '', args)
          args[[2]] <- get(selected_channel)
          args[[3]] <- paste(gating_method, get(selected_channel), sep="_")
        }
        
        ## guess the gate name if not supplied
        ## TODO: 2D gates?
        if (gate_name == '' && gating_method != "flowClust.2d") {
          if (args[["positive"]]) {
            suffix <- "+"
          } else {
            suffix <- "-"
          }
          gate_name <- paste0( get(selected_channel), suffix)
        }
        
        for (sample in samples) {
          
          args[[1]] <- quote(fs[[sample]])
          
          if (gate_both) {
            
            args[["positive"]] <- TRUE
            gate <- do.call(fn, args)
            flowWorkspace::add(gs[sample], gate, parent=parent, name=paste0( get(selected_channel), "+"))
            
            args[["positive"]] <- FALSE
            gate <- do.call(fn, args)
            flowWorkspace::add(gs[sample], gate, parent=parent, name=paste0( get(selected_channel), "-"))
            
          } else {
            gate <- do.call(fn, args)
            print(gate)
            flowWorkspace::add(gs[sample], gate, parent=parent, name=gate_name)
          }
          
        }
        
        recompute(gs[samples])
        
        ## update inputs to reflect changes
        allNodes <- unique( unlist( lapply(gs, function(x) {
          getNodes(x, isPath=TRUE)
        })))
        
        select <- grep( gate_name, allNodes, fixed=TRUE, value=TRUE )
        updateSelectInput(session, "parent_selection", choices=allNodes, selected=select)
        
        ## update the template data.frame to reflect the new gate
        gating_args <- args[4:length(args)]
        template <<- rbind(template, data.frame(
          alias=gate_name,
          pop=gate_name,
          parent=parent,
          dims=channel_1,
          gating_method=gating_method,
          gating_args=paste( names(gating_args), gating_args, sep="=", collapse="," ),
          collapseDataForGating='',
          groupBy='',
          preprocessing_method=if (grepl("flowClust", gating_method)) "prior_flowClust" else '',
          preprocessing_args='',
          samples=paste(samples, collapse=",")
        ) )
      }, error=function(e) {
        cat("ERROR: Could not apply gate!")
      })
        
    })
  })
  
  output$save_template <- downloadHandler(
    filename="GatingTemplate.csv",
    content=function(file) {
      write.table(file, row.names=FALSE, col.names=TRUE, sep=";", quote=FALSE)
    }
  )
  
  ## A GatingSet save observer
  output$save_gs <- downloadHandler(
    filename = "GatingSet.zip",
    content = function(file) {
      tempdir <- tempdir()
      dir <- file.path(tempdir, basename(tempfile()))
      zip_path <- file.path(tempdir, "gs.zip")
      save_gs(gs, path=dir)
      zip(zip_path, list.files(dir, full.names=TRUE))
      file.copy(zip_path, file)
      system( paste("rm -rf", dir) )
      file.remove(zip_path)
    },
    contentType="application/octet-stream"
  )
  
  output[["2d_densityplot"]] <- renderPlot({
    channel_1 <- input$channel_1
    channel_2 <- input$channel_2
    parent <- input$parent_selection
    samples <- getSampleSelection()
    
    ## only take the samples that have 'parent' as a node
    keep <- unlist( lapply( gs[samples], function(x) {
      return( parent %in% getNodes(x, isPath=TRUE) )
    } ) )
    
    samples <- samples[keep]
    if (length(samples) == 0) {
      grid.text("No information available at this node for the selected samples.")
      return( invisible(NULL) )
    }
    
    fs <- getData(gs[samples], parent)
    print( flowViz::xyplot( as.formula( paste(channel_2, "~", channel_1) ), fs ) )
    
  })
  
  output[["2d_densityplot_diag"]] <- renderUI({
    
    samples <- getSampleSelection()
    parent <- input$parent_selection
    missing_samples <- NULL
    
    for (sample in samples) {
      if (!(parent %in% getNodes(gs[[sample]], isPath=TRUE))) {
        missing_samples <- c(missing_samples, sample)
      }
    }
    
    if (!is.null(missing_samples)) {
      return( tags$div( class="gate-plot-diagnostics",
        h3("Samples Ungated at this Node"),
        p( paste(missing_samples, collapse=", ") )
      ) )
    } else {
      return( tags$div( class="gate-plot-diagnostics" ) )
    }
    
  })
  
  output[["1d_densityplot"]] <- renderPlot({
    
    channel_1 <- input$channel_1
    channel_2 <- input$channel_2
    n <- input[["1d_densityplot_n"]]
    parent <- input$parent_selection
    samples <- getSampleSelection()
    
    ## only take the samples that have 'parent' as a node
    keep <- unlist( lapply( gs[samples], function(x) {
      return( parent %in% getNodes(x, isPath=TRUE) )
    } ) )
    
    samples <- samples[keep]
    if (length(samples) == 0) {
      grid.text("No information available at this node for the selected samples.")
      return( invisible(NULL) )
    }
    
    fs <- getData(gs[samples], parent)
    
    ## subset for faster plotting
    fs <- fsApply(fs, function(x) {
      if (nrow(x) < n) {
        return (x)
      } else {
        return (x[sample(1:nrow(x), n), ])
      }
    })
    
    form <- as.formula(paste("~", channel_1, "+", channel_2))
    print( densityplot(form, fs) )
    
    return( invisible(NULL) )
    
  })
  
  output[["1d_densityplot_diag"]] <- renderUI({
    
    samples <- getSampleSelection()
    parent <- input$parent_selection
    missing_samples <- NULL
    
    for (sample in samples) {
      if (!(parent %in% getNodes(gs[[sample]], isPath=TRUE))) {
        missing_samples <- c(missing_samples, sample)
      }
    }
    
    if (!is.null(missing_samples)) {
      return( tags$div( class="gate-plot-diagnostics",
        h3("Samples Ungated at this Node"),
        p( paste(missing_samples, collapse=", ") )
      ) )
    } else {
      return( tags$div( class="gate-plot-diagnostics" ) )
    }
    
  })
  
  output$gate_plot <- renderPlot({
    samples <- getSampleSelection()
    parent <- input$parent_selection
    for (sample in samples) {
      if (!(parent %in% getNodes(gs[[sample]], isPath=TRUE))) {
        samples <- samples[ samples != sample ]
      }
    }
    
    if (length(samples)) {
      print( plotGate(gs[samples], parent, type="density") )
    }
    
    return( invisible(NULL) )
  })
  
  output$gate_plot_diag <- renderUI({
    
    samples <- input$sample_selection
    if (!length(samples)) {
      samples <- sampleNames(gs)
    }
    parent <- input$parent_selection
    missing_samples <- NULL
    
    for (sample in samples) {
      if (!(parent %in% getNodes(gs[[sample]], isPath=TRUE))) {
        missing_samples <- c(missing_samples, sample)
      }
    }
    
    if (!is.null(missing_samples)) {
      return( tags$div( class="gate-plot-diagnostics",
        h3("Samples with No Gates at this Node"),
        p( paste(missing_samples, collapse=", ") )
      ) )
    } else {
      return( tags$div( class="gate-plot-diagnostics" ) )
    }
    
  })
  

  output$gating_controls <- renderUI({
    gating_method <- input$gating_method
    fn <- get(gating_method, envir=asNamespace("openCyto"))
    fnArgs <- formals(fn)
    if (gating_method == "flowClust.2d") {
      fnArgs <- fnArgs[5:length(fnArgs)]
    } else {
      fnArgs <- fnArgs[4:length(fnArgs)]
    }
    fnArgs[sapply(fnArgs, is.null)] <- "NULL"
    fnArgs <- lapply(fnArgs, function(x) {
      if (is.call(x)) capture.output( print(x) )
      else as.character(x)
    })
    fnArgs <- fnArgs[ names(fnArgs) != "..." ]
    output <- vector("list", length(fnArgs))
    for (i in seq_along(fnArgs)) {
      arg <- fnArgs[[i]]
      nm <- names(fnArgs)[i]
      output[[i]] <- textInput(
        inputId=paste(gating_method, nm, sep="_"),
        label=nm,
        value=arg
      )
    }
    return(output)
  })
  
  output$debug <- renderPrint({
    
    R_send <- input$R_send
    if (input$R_send == 0) {
      return( invisible(NULL) )
    }
    
    isolate({
      code <- input$R_input
      result <- eval( parse( text=code ) )
      return(result)
    })
    
  })
  
  output$template <- renderTable({
    apply_gate <- input$apply_gate
    remove_gate <- input$remove_gate
    
    template
  })
  
  output$gating_hierarchy <- renderPlot({
    apply_gate <- input$apply_gate
    delete_gate <- input$delete_gate
    
    print( flowWorkspace::plot(gs[[1]]) )
  })
  
})
