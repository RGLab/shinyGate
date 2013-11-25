do_flowClust.2d <- function(env) {
  attach(env)
  ## the first three arguments are fixed and represent the flow frame,
  ## the channel name, and the filterId, for the 1d gates
  args$xChannel <- channel_1
  args$yChannel <- channel_2
  args$filterId <- paste(gating_method, paste(channel_1, channel_2, sep=":"), sep="_")
  args$lambda <- 1-1E-5
  
  fs <- getData(gs[samples], parent)
  
  ## default name if gate name unsupplied
  if (gate_alias == '') {
    gate_alias <- paste(channel_1, channel_2, sep=":")
  }
  
  ## if usePrior == "yes", we should generate a prior
  if (args[["usePrior"]] == "yes") {
    prior <- openCyto::prior_flowClust(fs, c(channel_1, channel_2), K=args[["K"]])
    args[["prior"]] <- prior
  }
  
  print(args)
  
  for (sample in samples) {
    
    args$fr <- substitute(fs[[sample]], list(sample=sample))
    gate <- do.call(fn, args)
    flowWorkspace::add(gs[[sample]], gate, parent=parent, name=gate_alias)
    
  }
  detach(env)
}
