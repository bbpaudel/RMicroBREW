
survivalRLS <- function(data_rls, numFrames, group_name, plot_it, save_file, color1){
  data_rls = data_rls
  data_rls = data_rls[data_rls$peak_height != 0, ]
  rls_data = data.frame()
  for(i in unique(data_rls$trap_id)){
    temp = subset(data_rls, data_rls$trap_id == i)
    count_data = table(temp$which_shape)
    count_data = data.frame(s1 = count_data[1], s2 = count_data[2], s3 = count_data[3])
    count_data$total = count_data$s1 + count_data$s2 + count_data$s3
    count_data$cell_entry = unique(temp$mother_cell_entry)
    count_data$death_frame = unique(temp$death_frame)
    count_data$live_status = ifelse(count_data$death_frame >= numFrames, 0, 1)
    count_data$trap_id = paste0(i)
    rls_data = rbind(rls_data, count_data)
  }
  rls_data$group = paste0(group_name)
  fit_rls = survfit(Surv(total, live_status) ~ group,
                 data = rls_data)
  print(fit_rls)
  if(plot_it == "yes"){
    print(ggsurvplot(fit_rls, data = rls_data,
               risk.table = F, surv.median.line = "none",xlim = c(0,50), break.time.by	= 10,
               combine = T, conf.int = T, pval = T, add.all = F, color = color1,
               xlab = "Generation", ylab = "Fraction Viable",
               risk.table.pos = 'in'))
  } else{
    print("figure not printed")
  }
  rls_data = rls_data[, c("trap_id", "cell_entry", "s1", "s2", "s3", "total", "death_frame", "live_status", "group")]

  outputdir = paste0(getwd(), "/RMicroBREW_results")
  if (!dir.exists(outputdir)) {dir.create(outputdir)}
  write.csv(rls_data, file = paste0(outputdir, "/", save_file, "RLS_count.csv"))
  return(rls_data)
}



