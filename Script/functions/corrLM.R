
corrLM <- function(data,x_col,y_col,x_lab="",title="",x_label,y_cor_label,y_formula_label){
  p <- ggplot(data, aes(x = get(x_col), y = get(y_col)))+
    geom_point(shape =3,color="gray60",size=0.7)+
    xlab(x_lab)+ylab(y_col)+labs(title=title)+
    theme_classic()+
    theme(title=element_text(size=9),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7))+
    geom_smooth(method ="lm",color ="black",fill ="lightblue3",size=0.6)+
    stat_cor(method ="spearman",cor.coef.name="rho",label.x =x_label, label.y=y_cor_label,size =2.5) +
    stat_poly_eq(
      aes(label = ..eq.label..),
      formula =y ~x,parse =TRUE, geom ="text",label.x =x_label,label.y =y_formula_label, hjust =0,size =2.5)
  return(p)
}





