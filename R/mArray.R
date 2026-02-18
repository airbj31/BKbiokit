#' generate dynamic methylation clock using glmnet package.
#'
#' @param x  tibble object which has n x m including SID, Age
#' @param y  Age column name
#' @param offset Gastrulation time. The value is used to construct aging value.
#' @param marker 1 dimentional character vector containing markers which passed the Quality Control
#' @param test test set array
#' @param useFemale T/F value
#' @param s         lambda.min or lambda.1se
#' @param alpha     1
#' @param weight    weight
#' @import glmnet
clock_gen <- function(x,y=Age,marker,test,offsets=0.05205479,SID=SID,useFemale=F,s=c("lambda.min","lambda.1se"),alpha=1, weight=1) {
  if(useFemale==T) {marker=c(marker,"Female")}
  yval <- logAge(x |> pull({{y}}),offsets)
  cvfit    <- cv.glmnet(x |> select(all_of(marker))  |> as.matrix(),yval,alpha=alpha,weight=weight)

  ## training result #1 lambda.min
  training.samples <- x |> pull({{SID}})

  cvresult<-list()
  for(i in s) {
    cvresult[[i]] <- predict(cvfit,newx=x |> select(all_of(marker))  |> as.matrix(),s=i)
    cvresult[[i]] <- tibble(SID=training.samples,Age=x |> pull({{y}}),pred.Age=expAge(cvresult[[i]][,1],offsets),method=i)
  }
  model_coef <- list()
  for(i in s) {
    model_coef[[i]] <- tibble(CGid=row.names(coef(cvfit,s=i)), beta=coef(cvfit,s=i)[,1]) |> dplyr::filter(beta>0)
  }
  MAE <- list()
  for(i in s) {
    MAE[[i]] <- abs(cvresult[[i]] |> pull({{y}}) - cvresult[[i]]$pred.Age) |> sum()
    MAE[[i]] <- MAE[[i]]/nrow(x)
  }
  model_stat <- list()
  for(i in s) {
    model_stat[[i]] <- glance(lm(cvresult[[i]] |> pull({{y}}) ~ pred.Age ,data=cvresult[[i]])) |> mutate(MAE=MAE[[i]],nCGid=nrow(model_coef[[i]]),s=i) |> select(s,everything())
  }
  model_stat <-do.call("rbind",model_stat)

  result <- list()
  for(i in s) {
    result[[i]] <- predict(cvfit,newx=as.matrix(test |> select(all_of(marker))),s=i)
    result[[i]] <- tibble(SID=row.names(result[[i]]),Age=test |> pull({{y}}),pred.Age=result[[i]][,1],method=i) |> mutate(pred.Age=expAge(pred.Age,offsets))
  }
  train.plot <- list()
  for(i in s) {
    tmp <- model_stat |> dplyr::filter(s==i) |> pull(adj.r.squared)
    train.plot[[i]] <- draw_train_result(cvresult[[i]],{{y}},scale="month",rsq=tmp,MAE=MAE[[i]])
  }
  return(list(fit=cvfit,
              model_coef=model_coef,
              model_stat=model_stat,
              train_result=cvresult,
              train_plot= train.plot,
              test_result=result
  ))
}

#'
RelAge <- function(Age,GestationT,maxLifeSpan) {
  RelAge <- (Age + GestationT)/(maxLifeSpan + GestationT)
  return(RelAge)
}

loglogAge <- function(x) {
  loglogAge = -log(-log(x))
  return(loglogAge)
}

DNAmAge <- function(x) {
  DNAmAge <- exp(-1 * ex(-1 *x)) * (maxLifespan + GestatT) - GestatT
  return(DNAmAge)
}

RelAdultAge <- function(Age,GestationT,ASM) {
  RelAdltAge <- (Age + GestationT)/(ASM + GestationT)
  return(RelAdltAge)
}

maxAgeParam <- function(c1,maxlifespan,gestationT,ASM) {
  m = c1 * (maxlifespan + gestationT)/(ASM + gestationT)
}


#' manipulate input age
#'
logAge <- function(x,offset) {
  return(log(x + offset))
}

#' manipulate output age
#'
expAge <- function(x,offset) {
  return(exp(x) - offset)

}

draw_train_result <- function(x,y,scale="year",rsq,MAE) {
  axis_tag <- paste0("(",scale,")")
  if(scale=="month") {scale <- 12} else {scale <- 1}
  train.plot <- ggplot(x, aes(x={{y}} * scale ,y=pred.Age * scale)) +
    geom_point(color="gray55",size=3) +
    geom_point(color="white",size=2.5) +
    geom_smooth(method="lm",alpha=0.5,color="gray33") +
    geom_richtext(data=tibble(Age=-Inf,pred.Age=Inf),
                  aes(label=paste0("adj-r<sup>2</sup> = ",round(rsq,digits=2),"<br>MAE = ",round(MAE,digits=3)),vjust=1,hjust=0,label.size=0)) +
    theme_bw() +  theme(legend.background=element_blank(),legend.position="none") +
    ylab(paste0("Predicted Age ",axis_tag)) +
    xlab(paste0("Age ",axis_tag))

  max_x <- max(x |> pull({{y}}))
  if(scale==1) {
    train.plot <- train.plot +
      scale_x_continuous(breaks=seq(0,max_x,by=0.25)) +
      scale_y_continuous(breaks=seq(0,max(x |> pull(pred.Age)),by=0.25))
  } else if(scale==12) {
    train.plot <- train.plot +
      scale_x_continuous(breaks=seq(0, max_x* 12,by=3)) +
      scale_y_continuous(breaks=seq(0,max(x |> pull(pred.Age) * 12),by=3))
  }
  train.plot <- train.plot + theme(panel.grid.minor=element_line(linetype="dashed"))
  return(train.plot)
}
