
# This file contains several functions used in plot for visulization

# utils::globalVariables(c("x", "value", "cv", "phi", "Phi", "pacf")) 




# plot the result of estimation 

plot_esti <- function(res_esti, ops, mp_type, title="", lower = -1.3, upper = 1.3, domain = 10){ # if fixt-->>2D plot fixt, if fixx  --->> fix time series plot,   nfix -->> 3D 
  
  n_esti = length(res_esti)
  
  if(ops == "fixt"){
    # 2D plot
    
    pos_map = c( "algebp", "logarip") 
    if(mp_type %in% pos_map){
      x_i = seq(0, domain, length.out = n_esti)
    }else{
      x_i = seq(-domain, domain, length.out = n_esti) 
    }
    
    df = data.frame(x=x_i, y = res_esti)
    
    theme_update(plot.title = element_text(hjust = 0.5))
    
    res = ggplot(df, aes_string(x="x", y= "y")) + geom_line(color = "#00AFBB")  + ggtitle(title) + ylim(lower, upper) +
      xlab("x") + ylab("m(t,x)")  + theme(plot.title = element_text(size=18, face="bold"),
                                                                             legend.text=element_text(size=24, face = "bold"),
                                                                             axis.text.x = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                             axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                             axis.title.x=element_text(size=22,face='bold'),
                                                                             axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                             legend.title = element_text(face = "bold"))
    # scale_colour_discrete(name  ="phi")
    return(res)
    
  } else if(ops == "fixx"){
    # 2D plot
    df = data.frame(t=seq(0, 1, length.out = n_esti), y = res_esti)
    
    theme_update(plot.title = element_text(hjust = 0.5))
    
    res = ggplot(df, aes_string(x="t", y="y")) + geom_line(color = "#00AFBB")  + ggtitle(title) + ylim(lower, upper) +
      xlab("t") + ylab("m(t,x)") + scale_colour_discrete(name  ="phi")+theme(plot.title = element_text(size=18, face="bold"),
                                                                             legend.text=element_text(size=24, face = "bold"),
                                                                             axis.text.x = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                             axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                             axis.title.x=element_text(size=22,face='bold'),
                                                                             axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                             legend.title = element_text(face = "bold"))
    return(res)
    
  } else if(ops == "nfix"){
    # 3D plot 
    fig <- plot_ly(z = ~res_esti)
    fig <- fig %>% add_surface()
    fig <- fig %>% layout(title = "Estimation function", scene = list(camera = list(eye = list(x = 2, y = 0.15, z = 1)),xaxis = list(title = 'x'),
                                                        yaxis = list(title = 't'),
                                                        zaxis = list(title = 'Este_function')))
    
    return(fig)
    
    
  } else{
    return(stop("Invalid option!"))
  }
}




# Cross validation plot,(with c, d and cross validation score)

fit.plot.cvm <- function(cv_m, title = ""){
  # library(plotly)
  df.cv = data.frame(c = cv_m[,1], b = as.factor(cv_m[,2]), cv = cv_m[,3])
  fig <- plot_ly(df.cv, x = ~c, y = ~b, z = ~cv,
                 marker = list(color = ~cv, colorscale = 'Viridis', showscale = TRUE),
                 text = ~paste('c:', c, '<br>b:', b, '<br>cv:', cv))
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(title = title, scene = list(camera = list(eye = list(x = -1.68, y = 1.68, z = 1.3)), xaxis = list(title = 'c'),
                                                    yaxis = list(title = 'b', tickvals = list(1,2)),
                                                    zaxis = list(title = 'cv')),
                        annotations = list(
                          x = 1.13,
                          y = 1.05,
                          text = 'cv',
                          xref = 'paper',
                          yref = 'paper',
                          showarrow = FALSE
                        ))
  return(fig)
  
}


# Plot SCR(including 2D)
plot_scr <- function(scr_df, ops, title = "", lower = -1.3, upper = 1.3){ # fixt-->>2D plot fixt, if fixx  --->> fix time series plot, scr_df is the result of SCR 
  if(ops == "fixt"){
    
    theme_update(plot.title = element_text(hjust = 0.5))
    

    res = ggplot(scr_df, aes_string(x="x", y="y", group = "order", colour = "order")) + geom_line() + ylim(lower, upper) + 
      scale_color_manual(values = c("#FF0000", "#3399FF", "#3399FF"), labels=c("EST", "LCR", "UCR"))  + ggtitle(title)  +
      xlab("x") + ylab("m(t,x)")  + guides(color = guide_legend(title = "")) +
                                          theme(plot.title = element_text(size=18, face="bold"),
                                          legend.text=element_text(size=22, face = "bold"),
                                          axis.text.x = element_text(face="bold", color="#993333",size=22, angle=0),
                                          axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                          axis.title.x=element_text(size=22,face='bold'),
                                          axis.title.y=element_text(angle=90, face='bold', size=22),
                                          legend.title = element_text(face = "bold"))
    # scale_colour_discrete(name  ="phi")
    return(res)
    
    
  } else if (ops == "fixx"){
    theme_update(plot.title = element_text(hjust = 0.5))
    res = ggplot(scr_df, aes_string(x = "t", y = "y", group = "order", colour = "order")) + geom_line() + ylim(lower, upper) + 
      scale_color_manual(values = c("#FF0000", "#3399FF", "#3399FF"), labels=c("EST", "LCR", "UCR"))  + ggtitle(title)  +
      xlab("t") + ylab("m(t,x)")  + guides(color = guide_legend(title = "")) +
                                          theme(plot.title = element_text(size=18, face="bold"),
                                          legend.text=element_text(size=22, face = "bold"),
                                          axis.text.x = element_text(face="bold", color="#993333",size=22, angle=0),
                                          axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                          axis.title.x=element_text(size=22,face='bold'),
                                          axis.title.y=element_text(angle=90, face='bold', size=22),
                                          legend.title = element_text(face = "bold"))
    # scale_colour_discrete(name  ="phi")
    return(res)
    
    
  } else{
    return(stop("Invalid option!"))
  }
  
}


# plot SCR

# plot time homogeniety
homot_plot <- function(sepera_df, title = "", lower = -1.3, upper = 1.3){
  theme_update(plot.title = element_text(hjust = 0.5))
  res = ggplot(sepera_df, aes_string(x="x", y="y", group = "order", colour = "order")) + geom_line() + ylim(lower, upper) +
    scale_color_manual(values = c("#FF0000", "#3399FF", "#3399FF"), labels=c("EST", "LCR", "UCR"))  + ggtitle(title)  +
    xlab("x") + ylab("m(t,x)")  + guides(color = guide_legend(title = "")) +
                                        theme(plot.title = element_text(size=18, face="bold"),
                                        legend.text=element_text(size=22, face = "bold"),
                                        axis.text.x = element_text(face="bold", color="#993333",size=22, angle=0),
                                        axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                        axis.title.x=element_text(size=22,face='bold'),
                                        axis.title.y=element_text(angle=90, face='bold', size=22),
                                        legend.title = element_text(face = "bold"))
  # scale_colour_discrete(name  ="phi")
  return(res)
  
}


# plot_sepera 
sep_plot <- function(sepera_df, ops, title = "", lower = -1.3, upper = 1.3){
  if(ops == "fixt"){
    
    theme_update(plot.title = element_text(hjust = 0.5))
    
    
    res = ggplot(sepera_df, aes_string(x="x", y="y", group = "order", colour = "order")) + geom_line() + ylim(lower, upper) + 
      scale_color_manual(values = c("#FF0000", "#3399FF", "#3399FF"), labels=c("EST", "LCR", "UCR"))  + ggtitle(title)  +
      xlab("x") + ylab("m(t,x)")  + guides(color = guide_legend(title = "")) + 
                                          theme(plot.title = element_text(size=18, face="bold"),
                                          legend.text=element_text(size=22, face = "bold"),
                                          axis.text.x = element_text(face="bold", color="#993333",size=22, angle=0),
                                          axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                          axis.title.x=element_text(size=22,face='bold'),
                                          axis.title.y=element_text(angle=90, face='bold', size=22),
                                          legend.title = element_text(face = "bold"))
    # scale_colour_discrete(name  ="phi")
    return(res)
    
  }else if(ops == "fixx"){
    theme_update(plot.title = element_text(hjust = 0.5))
    
    
    res = ggplot(sepera_df, aes_string(x="t", y="y", group = "order", colour = "order")) + geom_line() + ylim(lower, upper) +
      scale_color_manual(values = c("#FF0000", "#3399FF", "#3399FF"), labels=c("EST", "LCR", "UCR"))  + ggtitle(title)  +
      xlab("t") + ylab("m(t,x)")  + guides(color = guide_legend(title = "")) +
                                          theme(plot.title = element_text(size=18, face="bold"),
                                          legend.text=element_text(size=22, face = "bold"),
                                          axis.text.x = element_text(face="bold", color="#993333",size=22, angle=0),
                                          axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                          axis.title.x=element_text(size=22,face='bold'),
                                          axis.title.y=element_text(angle=90, face='bold', size=22),
                                          legend.title = element_text(face = "bold"))
    # scale_colour_discrete(name  ="phi")
    return(res)
    
  } else{
    return(stop("Invalid option!"))
  }
}

sep_plot_3d <- function(res, ops, title = title, lower = lower, upper = upper){
  
  t = seq(0,1, length.out = 300)
  x = seq(-10, 10, length.out = 300)
  plt_res = list()
  n_esti = dim(res[[2]][[1]])[1]
  k = length(res[[2]])
  for(i in 1:k){
    # 1 index for scr.df, 2 index for fix x point, 3 index for r 
    
    df_upper = res[[1]][[1]][[k]][res[[1]][[1]][[k]][,3] == "upper",c(1,2)]
    df_lower = res[[1]][[1]][[k]][res[[1]][[1]][[k]][,3] == "lower",c(1,2)]
    for(i in 2:n_esti){
      df_upper = rbind(df_upper, res[[1]][[i]][[k]][res[[1]][[i]][[k]][,3] == "upper",c(1,2)])
      df_lower = rbind(df_lower, res[[1]][[i]][[k]][res[[1]][[i]][[k]][,3] == "lower",c(1,2)])
    }
    df_upper$x = rep(x, each = n_esti)
    df_lower$x = rep(x, each = n_esti)
    # 1 index for scr.df, 2 index for fix x point, 3 index for r 
    
    aux_scr_upper = matrix(df_upper$y, ncol = n_esti, nrow = n_esti)
    aux_scr_lower = matrix(df_lower$y, ncol = n_esti, nrow = n_esti)
    estimate = res[[2]][[k]]
    
    fig <- plot_ly(showscale = FALSE) 
    fig <- fig %>% add_surface(x= ~t, y=~x, z = ~aux_scr_upper, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)")))
    fig <- fig %>% add_surface(x= ~t, y=~x, z = ~aux_scr_lower, colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)")))
    fig <- fig %>% add_surface(x= ~t, y=~x, z = ~estimate, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)")))
    
    fig <- fig %>% layout(title = title, scene = list(camera = list(eye = list(x = -1.68, y = 1.68, z = 1.3)), xaxis = list(title = 't'),
                                                      yaxis = list(title = 'x'),
                                                      zaxis = list(title = 'm_t_x')),
                          annotations = list(
                            x = 1.13,
                            y = 1.05,
                            text = 'm_t_x',
                            xref = 'paper',
                            yref = 'paper',
                            showarrow = FALSE
                          ))
    plt_res[[k]] = fig
  }
  
  return(plt_res)
}



# Code below should removed 