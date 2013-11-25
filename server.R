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
  GS <- gs <- load_gs("gs/LoveLab")
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
  
  source("common_functions.R")
  sapply( list.files("gating_wrappers", full.names=TRUE), source )
  
  ## Globals
  GS <- gs <- NULL
  FS <- fs <- NULL
  template <- NULL
  last_channel_index <- 1L
  
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
      
      ## update the parent selection
      allNodes <- unique( unlist( lapply(gs, function(x) getNodes(x, isPath=TRUE) ) ) )
      updateSelectInput(session, "parent_selection", choices=allNodes, selected=pp)
      
      ## auto-save
      isolate({
        if (input$autosave) {
          save_gs(gs, input$data_selection, overwrite=TRUE)
        }
      })
      
      ## update 2d overlay selections
      updateSelectInput(session, "2d_overlay", "Overlay", allNodes, "root")
      
    })
    
  })
  
  ## An observer for switching to different gs
  observe({
    debug("Trying to read and update select inputs.")
    file <- input$data_selection
    GS <<- gs <<- load_gs(file)
    FS <<- fs <<- getData(gs)
    channels <- sort(getData(gs)@colnames)
    nodes <- getNodes(gs[[1]], isPath=TRUE)
    updateSelectInput(session, "parent_selection", choices=nodes)
    updateSelectInput(session, "channel_1", choices=channels)
    updateSelectInput(session, "channel_2", choices=channels)
    updateSelectInput(session, "sample_selection", choices=sort(sampleNames(gs)))
    allNodes <- unique( unlist( lapply(gs, function(x) {
      getNodes(x, isPath=TRUE)
    })))
    updateSelectInput(session, "2d_overlay", "Overlay", allNodes, "root")
  })
  
  observe({
    
    selected_channel <- input$select_channels
    isolate({
      channel_1 <- input$channel_1
      if (selected_channel == "channel_1") {
        last_channel_index <<- 1L
      } else if (selected_channel == "channel_2") {
        last_channel_index <<- 2L
      }
    })
    
    
  })
  
  ## An observer for giving nice names to the radio selection.
  observe({
    
    channel_1 <- input$channel_1
    channel_2 <- input$channel_2
    
    isolate({
      selected_channel <- input$select_channels
      l <- setNames( list("channel_1", "channel_2"), c(channel_1, channel_2) )
      updateRadioButtons(session, 
        "select_channels", 
        "Select Channels to be used in Gating", 
        l,
        names(l)[[last_channel_index]]
      )
    })
    
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
      gate_alias <- input$gate_alias
      samples <- getSampleSelection()
      selected_channel <- input$select_channels
      gate_both <- input$gate_positive_and_negative
      
      ## only remove gates for samples that have gates at this node
      missing_samples <- NULL
      
      for (sample in samples) {
        if (!(parent %in% getNodes(gs[[sample]], isPath=TRUE))) {
          missing_samples <- c(missing_samples, sample)
        }
      }
      
      samples <- samples[ !(samples %in% missing_samples) ]
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
        txt <- input[[paste(gating_method, names(fnArgs)[i], sep="_")]]
        tmp <- eval(parse(text=txt))
        if (is.null(tmp)) {
          args[i] <- list(NULL)
        } else {
          args[[i]] <- tmp
        }
      }
      
      env <- environment()
      env$gs <- gs
      env$fs <- fs
      
      tryCatch({
        
        ## shiny likes to do terrible things to the execution stack,
        ## so we force these functions to evaluate in a given environment
        switch(gating_method,
          mindensity=do_mindensity(env),
          cytokine=do_cytokine(env),
          flowClust.1d=do_flowClust.1d(env),
          flowClust.2d=do_flowClust.2d(env),
          return("Could not select a gating method!")
        )
        
        rm(env)
        
        recompute(gs[samples])
        
        ## update inputs to reflect changes
        allNodes <- unique( unlist( lapply(gs, function(x) {
          getNodes(x, isPath=TRUE)
        })))
        
        select <- grep( gate_alias, allNodes, fixed=TRUE, value=TRUE )
        updateSelectInput(session, "parent_selection", choices=allNodes, selected=select)
        
        ## autosave?
        isolate({
          if (input$autosave) {
            save_gs(gs, input$data_selection, overwrite=TRUE)
          }
        })
        
        ## update the gates in the 2d gate select input
        updateSelectInput(session, "2d_overlay", "Overlay", allNodes, "root")
        
      }, error=function(e) {
        cat("ERROR: Could not apply gate!\n")
        cat("Message:\n--------\n\n")
        print(e)
      })
        
    })
  })
  
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
    overlay <- input[["2d_overlay"]]
    samples <- getSampleSelection()
    n <- input[["2d_densityplot_n"]]
    
    ## only take the samples that have 'overlay' as a node
    keep <- unlist( lapply( gs[samples], function(x) {
      return( overlay %in% getNodes(x, isPath=TRUE) )
    } ) )
    
    samples <- samples[keep]
    if (length(samples) == 0) {
      grid.text("No information available at this node for the selected samples.")
      return( invisible(NULL) )
    }
    
    #fs <- getData(gs[samples], parent)
    fs <- flowSet( lapply(gs[samples], function(x) {
      getData(x, parent)
    }))
    
    ## filter for faster plotting
    for (i in 1:length(fs)) {
      m <- nrow( fs[[i]] )
      if (m > n) {
        fs[[i]] <- fs[[i]][1:n, ]
      }
    }
    
    if (overlay != "root") {
      
      overlay_gate <- flowSet( lapply(gs[samples], function(x) {
        getData(x, overlay)
      }))
      
      print( flowViz::xyplot( as.formula( paste(channel_2, "~", channel_1) ), 
        fs, 
        smooth=FALSE,
        stats=TRUE,
        margin=TRUE,
        xbin=32,
        overlay=overlay_gate
      ) )
    } else {
      print( flowViz::xyplot( as.formula( paste(channel_2, "~", channel_1) ), 
        fs, 
        smooth=FALSE,
        stats=TRUE,
        margin=TRUE,
        xbin=32
      ) )
    }
    
    
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
    
    #fs <- getData(gs[samples], parent)
    fs <- flowSet( lapply( gs[samples], function(x) {
      getData(x, parent)
    }))
    
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
    
    ## should we plot a density plot?
    if (inherits( getGate( gs[[ samples[1] ]], parent ), "rectangleGate" )) {
      type <- "densityplot"
    } else {
      type <- "xyplot"
    }
    
    if (length(samples)) {
      print( plotGate(gs[samples], parent, type=type) )
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
  
  output$gating_hierarchy <- renderPlot({
    apply_gate <- input$apply_gate
    delete_gate <- input$delete_gate
    
    print( flowWorkspace::plot(gs[[1]]) )
  })
  
})
