do_flowClust.1d <- function(env) {
  attach(env)
  if (is.null(args[["K"]])) {
    warning("Changing 'K' from NULL to 2")
    args[["K"]] <- 2
  }
  
  fs <- getData(gs[samples], parent)
  
  ## the first three arguments are fixed and represent the flow frame,
  ## the channel name, and the filterId, for the 1d gates
  args <- c('', '', '', args)
  args[[2]] <- get(selected_channel)
  args[[3]] <- paste(gating_method, get(selected_channel), sep="_")
  
  
  ## guess the gate name if not supplied
  if (gate_alias == '') {
    if (args[["positive"]]) {
      suffix <- "+"
    } else {
      suffix <- "-"
    }
    gate_alias <- paste0( get(selected_channel), suffix)
  }
  
  for (sample in samples) {
    
    args[[1]] <- quote(fs[[sample]])
    
    if (gate_both) {
      
      args[["positive"]] <- TRUE
      gate <- do.call(fn, args)
      flowWorkspace::add(gs[sample], gate, parent=parent, name=paste0( gsub("[-\\+]", "", gate_alias), "+" ))
      
      args[["positive"]] <- FALSE
      gate <- do.call(fn, args)
      flowWorkspace::add(gs[sample], gate, parent=parent, name=paste0( gsub("[-\\+]", "", gate_alias), "-" ))
      
    } else {
      gate <- do.call(fn, args)
      print(gate)
      flowWorkspace::add(gs[sample], gate, parent=parent, name=gate_alias)
    }
    
  }
  detach(env)
}