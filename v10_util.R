################################################################################
#This file contains the util fonctions for ATE estimation
#in ADRENALINE project
################################################################################

#Load all necessary libraries
library(tidyverse)
library(ggsci)
library(survival)
library(survminer)
library(cowplot)
library(grDevices)
library(ggpubr)
library(finalfit)
library(gridExtra)
library(stringr)
library(gt)
library(gtsummary)
library(survey)
library(tableone)
library(cobalt)
library(WeightIt)
library(RISCA)
library(cowplot)
library(pracma)
library(writexl)
library(RColorBrewer)
library(forestplot)

folder_txt = "" #TO CHANGE

################################################################################
#Prepare dataset df for given molecule (ATC level5)
################################################################################
prepare_df_atc5 <- function(df_tot,atc_keep){

	#Create a variable med_bin : binarychronic intake for the given medication
	df_tot$med_bin <- df_tot[,paste0(atc_keep,"_bin")]	
	print(summary(factor(df_tot$med_bin)))
	
	#Select only necessary variables
	df <- df_tot %>% 
		select(cletri, year_diag, depriv_index, depriv_index_quintile,
				screening_before_diag, age_3cl,
				gp_visit_4cl, gyneco_visit_4cl, age, AGE_CL_3_CL, starts_with("comedic"), 				
				matches("^[A-Z]{1}[0-9]{2}_bin"),
				matches(paste0("^",substr(atc_keep,1,3),"[A-Z]{1}_bin")),
				matches(paste0("^",substr(atc_keep,1,4),"[A-Z]{1}_bin")),
				matches(paste0("^",substr(atc_keep,1,5),"[0-9]{2}_bin")),
				starts_with("comor"), pnuicc_2cl, subtype, subtype4, breast_surgery_3cl, 
				axillary_surgery, surg_medical_structure, ct, rt, ht, ht_type_3cl, 
				antiher2, primary_ttt_5cl, status_vital,status_vital_txt, delay_os, charlson_indx_age, med_bin) %>%			
		select(-paste0(substr(atc_keep,1,3),"_bin")) %>% #We suppress level2 ATC variable	
		select(-paste0(substr(atc_keep,1,4),"_bin")) %>% #We suppress level3 ATC variable
		select(-paste0(substr(atc_keep,1,5),"_bin")) %>% #We suppress level4 ATC variable
		select(-paste0(atc_keep,"_bin")) #We suppress level5 ATC variable
		
	#Final variable preprocessing (replace NA, relevel factors)
	df <- df %>% 
		mutate(depriv_index = ifelse(depriv_index == "No",NA, depriv_index)) %>%
		mutate(depriv_index_quintile = ifelse(depriv_index_quintile == "No",NA, depriv_index_quintile)) %>%
		mutate(comedic_count = case_when(
			med_bin == "Yes" ~ as.numeric(comedic_count) -1,
			TRUE ~ as.numeric(comedic_count)
		)) %>%
		mutate(comedic_count_cl = cut(as.numeric(comedic_count),
							c(-1,0,5,11,Inf), 
							right = TRUE,
							labels = c("0","1-5","6-11","12+"))) %>%
		mutate(var = factor(med_bin, levels = c("Yes","No")),
				delay_os =  as.numeric(delay_os),
				status_vital = as.numeric(status_vital), 
				med_numeric = ifelse(med_bin == "Yes",1,0),
				med_bin = factor(med_bin),
				depriv_index = as.numeric(depriv_index), 
				age = as.numeric(age),
				gp_visit_4cl = factor(gp_visit_4cl, levels = c("0","1-5","6-11","12+")),
				gyneco_visit_4cl = factor(gyneco_visit_4cl, levels = c("0","1","2-3","4+")), 
				richness_index = -as.numeric(depriv_index),
				age_3cl = factor(age_3cl, levels = c("<50","50-75",">75")),
				age_cl_3_cl = cut(age, c(0,50,65,Inf), 
										right = FALSE, 
										labels = c("[0-50[","[50-65[","65+")))
										
	return(df)									

}

################################################################################
  # First row of plots : 
  # - bar plot for comedication intake yes versus no
  # - violin plot for age by comedication intake
  #  -violin plot for richness index by comedication intake
  # -percentage of patients with at least one comorbid condition by comedication intake
  # -percentage of obsese patients with  condition by comedication intake
  
  # Plot are saved individually in the folder 'descriptive'
  # And merged all together in the folder 'plot_merged'
################################################################################

plot_row1_boxplot <- function(id,atc,atc_name,atc_color,
											atc_level = "atc2",
											estimate = "ATE", label = c("A","B")){

	#Dataset preprocessing
	id_comor <- id %>% 
		select(cletri,var) %>%
		mutate(n_total = length(unique(cletri))) %>%
		group_by(var) %>%
		summarise(count = n(), percent = 100*count/n_total) %>%
		unique() %>%
		mutate(label = paste0("n = ",format(count,big.mark  = " "), "\n(", round(percent,1),"%)"))
	
	###################################
	#Bar plot for comedication intake yes versus no
	###################################
	#Create plot
	bar_use <- ggplot(data = id_comor, aes(x = var, y = percent, fill = var)) +
		geom_bar(stat = "identity") + 
		geom_text(aes(label = label),vjust = -0.5, size = 3) +  
		theme_classic() + 
		theme(legend.position = 'none') + 
		scale_y_continuous(limits = c(0,1.2*max(id_comor$percent)))+
		labs(subtitle = paste0(atc," consumption")) +
		ggtitle(atc) + 
		scale_fill_manual(values = c(atc_color,"#868686")) + 
		labs(x = atc_name, y = "Percentage of patients") + 
		theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5))) +
		theme(plot.title = element_text(size = rel(1.2), color = "black"))
	#Save plot individually
	png(paste0(folder_txt,atc_level,"/",estimate,"/descriptive/",atc,"_page1A.png"),
		height = 2000, width = 2000, units = "px", res = 300,
		)
		print(bar_use)
	dev.off()
	
	###################################
	#Violin plot for age by comedication intake
	###################################
	#Create plot
	subtitle_str = paste0("Age at diagnostic\nn = ", format(nrow(id), big.mark =  " "), " (100%)")
		#Box_plot_comed is an internal function
	box_age <- box_plot_comed(id,var_cont = "age",subtitle_str  = subtitle_str) 
	#Save plot
	png(paste0(folder_txt,atc_level,"/",estimate,"/descriptive/",atc,"_page1B.png"),
		height = 2000, width = 2000, units = "px", res = 300,
		)
		print(box_age)
	dev.off()
	
	########################################
	#Violin plot for richness index by comedication intake
	########################################	
	#Create plot
	n_depriv <- length(which(!is.na(id$richness_index)))
	subtitle_str = paste0("Richness index\nn = ", format(n_depriv, big.mark = " "), " (", round(100*n_depriv/nrow(id),1),"%)")
		#Box_plot_comed is an internal function
	box_age_depriv <- box_plot_comed(id,var_cont = "richness_index",subtitle_str  = subtitle_str, title = FALSE) 
	#Save plot	
	png(paste0(folder_txt,atc_level,"/",estimate,"/descriptive/",atc,"_page1C.png"),
		height = 2000, width = 2000, units = "px", res = 300,
		)
		print(box_age_depriv)
	dev.off()

	#Plot for any comorbidity and obesity 
	#Bar plot comorbidites (sans obesity et sans any comorbidity)
	id_comor <- id %>% 
		select(cletri,comorbidity_bin,comor2_Obesity, var) %>%
		pivot_longer(cols = starts_with("comor"), names_to = "comor", names_prefix = "", values_to = "bin") %>%
		mutate(n_total_pop = length(unique(cletri))) %>% 
		group_by(var) %>%
		mutate(n_total = length(unique(cletri))) %>%
		mutate(comor = factor(comor, levels = c("comorbidity_bin","comor2_Obesity"), 
									 labels = c("Any comorbidity","Proxy for obesity")		 
							)) %>%
		group_by(comor) %>%
		mutate(count_comor = length(unique(cletri[bin == "Yes"])), 
			   percent_comor = 100*count_comor/n_total_pop,
			   facet_label = paste0(comor,"\nn = ",format(count_comor,big.mark  = " "), " (", round(percent_comor,1),"%)")) %>%
		group_by(comor,var,facet_label) %>%
		summarise(count = length(unique(cletri[bin == "Yes"])), percent = 100*count/n_total) %>%
		unique() %>%
		mutate(label = paste0("n = ",format(count,big.mark  = " "), "\n(", round(percent,1),"%)"))
	
	bar_comor <- ggplot(data = id_comor, aes(x = var, y = percent, fill = var)) +
		geom_bar(stat = "identity") + 
		geom_text(aes(label = label),vjust = -0.5, size = 3) +  
		theme_classic() + 
		theme(legend.position = 'none') + 
		scale_y_continuous(limits = c(0,1.2*max(id_comor$percent)))+
		scale_fill_manual(values = c(atc_color,"#868686")) + 
		labs(x = atc_name, y = "Percentage of patients") +
		facet_wrap(.  ~ facet_label, nrow = 1)  +
		ggtitle(paste0("Distribution of comorbidities and proxy for morbid obesity\nw.r.t. ", atc_name, " use" )) +
		theme(plot.title = element_text(size = rel(1.2), color =  "black"))
	
	png(paste0(folder_txt,atc_level,"/",estimate,"/descriptive/",atc,"_page1D.png"),
		height = 2000, width = 4000, units = "px", res = 300,
		)
		print(bar_comor)
	dev.off()

	#Plot grid row 1 
	row1 <- plot_grid(bar_use,box_age,box_age_depriv, bar_comor, labels = c(label[1],'','',label[2]), nrow = 1, rel_widths = c(1,1,1,2))
	return(row1)
}

################################################################################
#Internal function
#Box plot for continusous variable (var_cont) by comedication use
################################################################################
box_plot_comed <- function(df, var_cont,subtitle_str, title = TRUE){	
	
	df[,"var_cont"] <- df[, var_cont]
	df$var_cont <- as.numeric(as.character(df$var_cont))
	
	title_str <- paste0(
				atc_name,
				"\n",
				"n = ",
				format(sum(df$var == "Yes"),big.mark = " "),
				" / ", 
				format(nrow(df),big.mark =" "),
				" (",
				round(100*sum(df$var == "Yes")/nrow(df),1) ,"%)"
			)
	
	box_age <- ggplot(data = df, aes(x = var, y = var_cont, fill = var)) +
		geom_violin() + 
		scale_fill_manual(values = c(atc_color,"#868686")) +
		xlab("")+
		ylab("") + 
		labs(subtitle = subtitle_str) +
		scale_y_continuous(expand = c(0.1,0.5)) + 
		ggtitle(title_str) + #peut etre mettre en blanc pour meilleur alignement
		stat_compare_means(method = "t.test") +
		theme_classic() +
		theme(legend.position = 'none') + 
		theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5))) +
		theme(plot.title = element_text(size = rel(1.2), color = ifelse(title, "black","white")))
		
	return(box_age)
			
}

################################################################################
#Descriptive plots row 2 to 5
#Bar plots for all comorbid condition by medication intake except obesity (already in row 1)
################################################################################

plot_row2345_boxplot <- function(id,atc_name,atc_color,
												  atc_level  ="atc2",
												  estimate = "ATE", label = c("C")){

	#Final data preprocesing for the plot
	id_comor <- id %>% 
		select(cletri,starts_with("comor2"), var) %>%
		select(- comor2_Obesity) %>% 
		pivot_longer(cols = starts_with("comor"),
							names_to = "comor",
							names_prefix = "comor2_",
							values_to = "bin") %>%
		mutate(n_total_pop = length(unique(cletri))) %>% 
		group_by(var) %>%
		mutate(n_total = length(unique(cletri))) %>%
		filter(!is.na(comor)) %>% 
		mutate(comor = gsub("_" ," ", comor)) %>%
		group_by(comor) %>%
		mutate(count_comor = length(unique(cletri[bin == "Yes"])), 
			   percent_comor = 100*count_comor/n_total_pop,
			   facet_label = paste0(comor,"\nn = ",format(count_comor,big.mark  = " "), " (", round(percent_comor,1),"%)")) %>%
		group_by(comor,var,facet_label) %>%
		summarise(count = length(unique(cletri[bin == "Yes"])), percent = 100*count/n_total) %>%
		unique() %>%
		mutate(label = paste0("n = ",format(count,big.mark  = " "), "\n(", round(percent,1),"%)"))	
										          
	#Order comorbidites by categories
	levels_label_facet <- levels(as.factor(id_comor$facet_label))
	levels_reorder <- c(levels_label_facet[grepl("Hypertension",levels_label_facet)],
						levels_label_facet[grepl("Congestive",levels_label_facet)],
						levels_label_facet[grepl("Ischemic",levels_label_facet)],
						levels_label_facet[grepl("Arrhythmia",levels_label_facet)],
						levels_label_facet[grepl("Valvular",levels_label_facet)],
						levels_label_facet[grepl("Cerebrovascular",levels_label_facet)],
						levels_label_facet[grepl("embolism",levels_label_facet)],
						levels_label_facet[grepl("Myocardiopathy",levels_label_facet)],
						levels_label_facet[grepl("Haemorrhagic",levels_label_facet)],
						levels_label_facet[grepl("Periph. vasc.",levels_label_facet)],							
						levels_label_facet[grepl("Diabete",levels_label_facet)],						 
						levels_label_facet[grepl("Thyroid",levels_label_facet)],														
						levels_label_facet[grepl("kidney",levels_label_facet)],
						levels_label_facet[grepl("Liver",levels_label_facet)],
						levels_label_facet[grepl("Pancreatitis",levels_label_facet)],
						levels_label_facet[grepl("Peptic",levels_label_facet)],
						levels_label_facet[grepl("Inflammatory",levels_label_facet)],
						levels_label_facet[grepl("HIV",levels_label_facet)],												
						levels_label_facet[grepl("Epilepsy",levels_label_facet)],					
						levels_label_facet[grepl("Dementia",levels_label_facet)], 
						levels_label_facet[grepl("Parkinson",levels_label_facet)], 
						levels_label_facet[grepl("sclerosis",levels_label_facet)],					 				
						levels_label_facet[grepl("Psychiatric",levels_label_facet)],
						levels_label_facet[grepl("Cognitive",levels_label_facet)],
						levels_label_facet[grepl("Down",levels_label_facet)],					
						levels_label_facet[grepl("Frailty",levels_label_facet)], 			
						levels_label_facet[grepl("Smoking",levels_label_facet)], 					
						levels_label_facet[grepl("Alcohol",levels_label_facet)],
						levels_label_facet[grepl("substance",levels_label_facet)],						
						levels_label_facet[grepl("Respiratory",levels_label_facet)],				
						levels_label_facet[grepl("Rheumatologic",levels_label_facet)],
						levels_label_facet[grepl("Connective",levels_label_facet)],
						levels_label_facet[grepl("Muscle disorders",levels_label_facet)],					
						levels_label_facet[grepl("Metabolic",levels_label_facet)],
						levels_label_facet[grepl("Transplan",levels_label_facet)],
						levels_label_facet[grepl("paralysis",levels_label_facet)]
						)

	id_comor_desc <- id_comor %>% 
		group_by(facet_label) %>% 
		summarise(count = sum(count)) %>%
		arrange(desc(count)) %>%
		ungroup()
	id_comor <- id_comor %>%
		mutate(facet_label = factor(facet_label, levels= unique(id_comor_desc$facet_label))) 
	
	#All bar plots 
	bar_comor <- ggplot(data = id_comor, aes(x = var, y = percent, fill = var)) +
		geom_bar(stat = "identity") + 
		geom_text(aes(label = label),hjust = -0.3, size = 4) +  
		theme_classic() + 
		theme(legend.position = 'none') + 
		scale_y_continuous(limits = c(0,1.7*max(id_comor$percent)))+
		scale_fill_manual(values = c(atc_color,"#868686")) + 
		labs(x = atc_name, y = "Percentage of patients") +
		coord_flip() + 
		facet_wrap(.  ~ facet_label, nrow = 10)  +
		ggtitle(paste0("Distribution of comorbidities w.r.t. ", atc_name, " use" )) + 
		theme(plot.title = element_text(size = rel(1.2), color =  "black"))
	
	#Save individual plot	
	png(paste0(folder_txt,atc_level,"/",estimate,"/descriptive/",atc,"_page1E.png"),
		height = 5333, width = 6000, units = "px", res = 300,
		)
		print(bar_comor)
	dev.off()
	
	#Add label	
	row2345 <- plot_grid(bar_comor, labels = c(label[1]), nrow = 1)
 	
	return(row2345)

}

################################################################################
#Compute IPTW weights
################################################################################
compute_iptw_weights <- function(id, df_comor, col_comed) {

	#List of variables included in the model
	explanatory_simple = c("year_diag", "age" , "screening_before_diag",
			"depriv_index_quintile",
			"gp_visit_4cl","gyneco_visit_4cl", 
			"comedic_count_cl",
		    col_comed,
		    colnames(df_comor))
	explanatory_simple_txt <- paste0(explanatory_simple, collapse = " + ")
	
	#Simple logisitc regression with weighted outputs: compute weights.
	tweights <- table(id$med_numeric)
	id$weights <- ifelse(id$med_numeric == 0,(500)/tweights[1],(500)/tweights[2]) 
	p_logistic.fit_text <- paste0("glm( med_numeric ~ ", explanatory_simple_txt, ", data = id, family = binomial(link = 'logit'), weights= id$weights)")
	p_logistic.fit <- eval(parse(text = p_logistic.fit_text))
	
	#Propensity score computation from logistic models
	id$p_logistic_weights <- as.vector(predict(p_logistic.fit, type = "response"))
			
	#IPT weights computation for propensity scores
	#Weights are trimmed
	id <- id %>% mutate(
		freq_cel = sum(med_numeric == 1)/n(), 
		IPW_logistic_stabilized = ifelse(med_numeric == 0, (1-freq_cel) /(1 - p_logistic_weights), freq_cel / p_logistic_weights),
		quantile_95 = quantile(IPW_logistic_stabilized, 0.975),
		quantile_05 = quantile(IPW_logistic_stabilized, 0.025),
		IPW_logistic_stab_trunc_right  = ifelse(IPW_logistic_stabilized > quantile_95, quantile_95, IPW_logistic_stabilized ))
	
	return(id)
}

################################################################################
#Compute adjustment quality
################################################################################
compute_adjustment_quality <- function(id, df_comor, col_comed, atc, atc_level="atc2", estimate = "ATE") {

	#Dataframe with pretty colnames
	pretty_colnames <- readxl::read_xlsx(paste0(folder,"ref/pretty_colnames.xlsx"))
	
	#######################################
	#First set of variables : socio-demographic variables
	#######################################
	
	#Prepare list of variables and pretty colnames
	explanatory1 = c( "age", "screening_before_diag",
					"depriv_index_quintile",
					"gp_visit_4cl","gyneco_visit_4cl", 
					"comedic_count_cl")				
	covs <- id[, explanatory1]
	colnames_pretty <- data.frame(old_colnames = explanatory1) %>% left_join(pretty_colnames)
	colnames(covs) <- colnames_pretty$pretty_colnames
	
	#Compute SMD
	W <- as.weightit(weights = id$IPW_logistic_stab_trunc_right, treat = id$med_numeric, c = covs)
	tab1 <- bal.tab(W,un = TRUE, m.threshold = 0.1, v.threshold =2)
	love_plot1 <- do.call(love.plot, args = list(x = tab1,colors =c("#1D3557","#E63946")))
	tab1_max <- tab1$Max.Imbalance.mean.diffs
	
	#Save love plot individually
	png(paste0(folder_txt,atc_level,"/", estimate, "/adjustment/",atc,"_page2B.png"),
		height = 2000, width = 2000, units = "px", res = 300,
		)
		print(love_plot1)
	dev.off()
	
	#######################################
	#Second set of variables : comorbid conditions part 1
	#######################################
	
	#Prepare list of variables and pretty colnames
	explanatory2 = c(colnames(df_comor)[1:25])					
	covs <- id[, explanatory2]
	colnames_pretty <- data.frame(old_colnames = explanatory2) %>% left_join(pretty_colnames)
	colnames(covs) <- colnames_pretty$pretty_colnames

	#Compute SMD
	W <- as.weightit(weights = id$IPW_logistic_stab_trunc_right, treat = id$med_numeric, covs = covs)
	tab2 <- bal.tab(W,un = TRUE, m.threshold = 0.1, v.threshold =2)
	love_plot2 <- do.call(love.plot, args = list(x = tab2,colors =c("#1D3557","#E63946")))
	tab2_max <- tab2$Max.Imbalance.mean.diffs
	
	#Save love plot individually
	png(paste0(folder_txt,atc_level,"/",estimate,"/adjustment/",atc,"_page2C.png"),
		height = 2000, width = 2000, units = "px", res = 300,
		)
		print(love_plot2)
	dev.off()
	
	#######################################
	#Third set of variables : comorbid conditions part 2
	#######################################
	
	#Prepare list of variables and pretty colnames
	explanatory3 = c(colnames(df_comor)[26:length(colnames(df_comor))])					
	covs <- id[, explanatory3]
	colnames_pretty <- data.frame(old_colnames = explanatory3) %>% left_join(pretty_colnames)
	colnames(covs) <- colnames_pretty$pretty_colnames
	
	#Compute SMD
	W <- as.weightit(weights = id$IPW_logistic_stab_trunc_right, treat = id$med_numeric, covs = covs)
	tab3 <- bal.tab(W,un = TRUE, m.threshold = 0.1)
	love_plot3 <- do.call(love.plot, args = list(x = tab3,colors =c("#1D3557","#E63946")))
	tab3_max <- tab3$Max.Imbalance.mean.diffs
	
	#Save love plot individually
	png(paste0(folder_txt,atc_level,"/",estimate,"/adjustment/",atc,"_page2D.png"),
		height = 2000, width = 2000, units = "px", res = 300,
		)
		print(love_plot3)
	dev.off()
	
	#######################################
	#Fourth set of variables : other comedications part 1
	#######################################
	
	#Prepare list of variables and pretty colnames
	n_comed = length(col_comed)
	break_comed = as.integer(n_comed/2)
	explanatory4 = col_comed[1:break_comed]			
	covs <- id[, explanatory4]
	colnames_pretty <- data.frame(old_colnames = explanatory4) %>% left_join(pretty_colnames)
	colnames(covs) <- colnames_pretty$pretty_colnames
	
	#Compute SMD
	W <- as.weightit(weights = id$IPW_logistic_stab_trunc_right, treat = id$med_numeric, covs = covs)
	tab4 <- bal.tab(W,un = TRUE, m.threshold = 0.1, v.threshold =2)
	love_plot4 <- do.call(love.plot, args = list(x = tab4,colors =c("#1D3557","#E63946")))
	tab4_max <- tab4$Max.Imbalance.mean.diffs
	
	#Save love plot individually
	png(paste0(folder_txt,atc_level,"/",estimate,"/adjustment/",atc,"_page2E.png"),
		height = 2000, width = 2000, units = "px", res = 300,
		)
		print(love_plot4)
	dev.off()
	
	#######################################
	#Fifth set of variables : other comedications part 2
	#######################################
	
	#Prepare list of variables and pretty colnames	
	explanatory5 = col_comed[(break_comed+1):n_comed]				
	covs <- id[, explanatory5]
	colnames_pretty <- data.frame(old_colnames = explanatory5) %>%
		left_join(pretty_colnames) %>%
		mutate(pretty_colnames = ifelse(
				is.na(pretty_colnames),
				gsub("_bin","",old_colnames),
				pretty_colnames
			)
		)
	colnames(covs) <- colnames_pretty$pretty_colnames
	
	#Compute SMD
	W <- as.weightit(weights = id$IPW_logistic_stab_trunc_right, treat = id$med_numeric, covs = covs)
	tab5 <- bal.tab(W,un = TRUE, m.threshold = 0.1, v.threshold =2)
	love_plot5<- do.call(love.plot, args = list(x = tab5,colors =c("#1D3557","#E63946")))
	tab5_max <- tab5$Max.Imbalance.mean.diffs
	
	#Save love plot individually
	png(paste0(folder_txt,atc_level,"/",estimate,"/adjustment/",atc,"_page2F.png"),
		height = 2000, width = 2000, units = "px", res = 300,
		)
		print(love_plot5)
	dev.off()

	#######################################
	#Bind all data frames and summarise
	#######################################
	tab_complete <- rbind(tab1_max,tab2_max,tab3_max,tab4_max,tab5_max)
	tab_max <- tab_complete %>% as.data.frame() %>% mutate(Diff.Adj_abs = abs(Diff.Adj)) %>%  filter(Diff.Adj_abs== max(Diff.Adj_abs))
	tab_complete <- rbind(as.data.frame(tab1$Balance %>% select(rmst_unadjusted = Diff.Un, rmst_adjusted = Diff.Adj)),
					as.data.frame(tab2$Balance %>% select(rmst_unadjusted = Diff.Un, rmst_adjusted = Diff.Adj)),
					as.data.frame(tab3$Balance %>% select(rmst_unadjusted = Diff.Un, rmst_adjusted = Diff.Adj)),
					as.data.frame(tab4$Balance %>% select(rmst_unadjusted = Diff.Un, rmst_adjusted = Diff.Adj)),
					as.data.frame(tab5$Balance %>% select(rmst_unadjusted = Diff.Un, rmst_adjusted = Diff.Adj)))
					
	return(list("love_plot1" = love_plot1,
				"love_plot2" = love_plot2,
				"love_plot3" = love_plot3,
				"love_plot4" = love_plot4,
				"love_plot5" = love_plot5,
				"tab_complete" = tab_complete,
				"tab_max"    =  tab_max))
}

################################################################################
#Compute density plots for propensity score
################################################################################
plot_ps <- function(id,name_atc_tokeep,atc_color,atc, atc_level=  "atc2", estimate = "ATE", print = TRUE) {

	#Density plot
	plot_weights <- ggplot(data =id, aes(x = p_logistic_weights , color = med_bin)) + 
			geom_density(size = 1) + 
			scale_color_manual(values = c("Yes"=atc_color,"No"="#868686")) + 
			theme_classic() + 
			xlab("Propensity score") +
			ylab("Density of probability") +
			labs(title =  paste0( "Distribution of propensity score"),
				color = name_atc_tokeep) +
			theme(legend.position = "bottom") 
	
	#Save as individual plot
	if(print)	{
		png(paste0(folder_txt,atc_level,"/",estimate,"/adjustment/",atc,"_page2A.png"),
			height = 2000, width = 2000, units = "px", res = 300,
			)
			print(plot_weights)
	dev.off()
}
	return(plot_weights)
}


################################################################################
#Kaplan-Meier survival curves WP, chemotherapy, no chemotherapy
################################################################################

plot_survival_wp_ct_no_ct <- function(id,atc,atc_name,atc_color, data_ct, data_no_ct, label1 = 'D', label2 = 'E',atc_level ="atc2",estimate  ="ATE",
											ylab_txt = "Overall survival probability"){
		
	print(dim(id))
	text_all <- case_when(
		estimate == "ATT"  ~"Treated population",
		estimate == "matched" ~ "Matched population",
		TRUE ~ "Whole population"
	)
	label_facet = paste0(text_all, ", n = ",format(nrow(id), big.mark = " ")," (100%)")
	#WP to have associated risk table 
	survA <- survfit(Surv(delay_os, status_vital) ~ var , data = id)
	plotA <- ggsurvplot(fit = survA,
						data = id,
						risk.table = TRUE,
						pval = TRUE,
						ggtheme = theme_bw(), 
						palette = c(atc_color,"#868686"),
						risk.table.col = "strata",
						risk.table.title = "Number at risk",
						legend = "none",
						legend.title = "",
						legend.labs = c("Yes","No"),
						censor = F,
						ylab = ylab_txt,
						risk.table.fontsize = 3.5)
	plotA$table

	#Survfit KM weighted
	survAKM_weighted <- ipw.survival(times = id$delay_os, failures = id$status_vital,
									variable = id$med_numeric, weights = id$IPW_logistic_stab_trunc_right)
								
	pA_value <- ipw.log.rank(times = id$delay_os, failures = id$status_vital,
									variable = id$med_numeric, weights = id$IPW_logistic_stab_trunc_right)
	print(pA_value)
	
	survAKM_weighted_df <- survAKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)	
		
	survAKM_weighted_df$strata <- factor(survAKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
	p_val_toprint <- ifelse(pA_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pA_value$p.value,3))))
	print(p_val_toprint)
	
	survA_plot <- ggsurvplot(fit = survAKM_weighted_df,
							palette = c(atc_color,"#868686"),
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
	survA_plot <- survA_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("Mortality by ",ifelse(nchar(atc_name) > 20, atc,atc_name)))),
							   atop(.(paste0("ATC Code : ", atc)))))) +
				labs(subtitle = label_facet) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold"))
	survA_plot   

#Part 3 : plot B : survival - CT Yes/No

		#Test interaction CT vs no CT
		mod1 <- coxph(Surv(delay_os, status_vital) ~  ct + var + ct*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
		mod2 <- coxph(Surv(delay_os, status_vital) ~ ct + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
		tmp_anova <- anova(mod1,mod2,test = "Chisq") 
		label_interaction <- paste0("Interaction ",atc,"-CT p-value : ", 
								ifelse(round(tmp_anova[,"P(>|Chi|)"][2],3) == 0, "<0.001", round(tmp_anova[,"P(>|Chi|)"][2],3))) 
		print(label_interaction)
  
 			#plots
		label_ct <- paste0("Chemotherapy, n = ", format(nrow(data_ct),big.mark = " "), " (", round(100*nrow(data_ct)/nrow(id),1),"%)")
		#Surv raw for risk table
		survB_ct <- survfit(Surv(delay_os, status_vital) ~ var, data = data_ct)
		pval_toprint_ct = surv_pvalue(survB_ct)$pval.txt
		plotB_ct <- ggsurvplot(fit = survB_ct, 
								data = data_ct,
								risk.table = TRUE, 
								pval = pval_toprint_ct, 
								ggtheme = theme_bw(), 
								palette = c(atc_color,"#868686"),
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = c("Yes","No"),
								censor = FALSE,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5) 
					
	#Survfit KM weighted
	survB_ctKM_weighted <- ipw.survival(times = data_ct$delay_os, failures = data_ct$status_vital,
									variable = data_ct$med_numeric, weights = data_ct$IPW_logistic_stab_trunc_right)
								
	pB_ct_value <- ipw.log.rank(times = data_ct$delay_os, failures = data_ct$status_vital,
									variable = data_ct$med_numeric, weights = data_ct$IPW_logistic_stab_trunc_right)
	print(pB_ct_value)
	
	survB_ctKM_weighted_df <- survB_ctKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
	head(survB_ctKM_weighted_df)	
		
	survB_ctKM_weighted_df$strata <- factor(survB_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
	p_val_toprint <- ifelse(pB_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pB_ct_value$p.value,3))))
	print(p_val_toprint)
	
	survB_ct_plot <- ggsurvplot(fit = survB_ctKM_weighted_df,
							palette = c(atc_color,"#868686"),
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
	survB_ct_plot <- survB_ct_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("Mortality by ",ifelse(nchar(atc_name) > 20, atc,atc_name), " and CT status"))),
							   atop(.(label_interaction))))) +
				labs(subtitle = label_ct) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold"))
	survB_ct_plot  
	
	#No CT
	label_no_ct <- paste0("No chemotherapy, n = ", format(nrow(data_no_ct),big.mark = " "), " (", round(100*nrow(data_no_ct)/nrow(id),1),"%)")
	#Surv raw for risk table
		survB_no_ct <- survfit(Surv(delay_os, status_vital) ~ var, data = data_no_ct)
		pval_toprint_no_ct = surv_pvalue(survB_no_ct)$pval.txt
		plotB_no_ct <- ggsurvplot(fit = survB_no_ct, 
								data = data_ct,
								risk.table = TRUE, 
								pval = pval_toprint_no_ct, 
								ggtheme = theme_bw(), 
								palette = c(atc_color,"#868686"),
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = c("Yes","No"),
								ylab = ylab_txt,
								censor = FALSE,
								risk.table.fontsize = 3.5) 
	
			#Survfit KM weighted
	survB_no_ctKM_weighted <- ipw.survival(times = data_no_ct$delay_os, failures = data_no_ct$status_vital,
									variable = data_no_ct$med_numeric, weights = data_no_ct$IPW_logistic_stab_trunc_right)
								
	pB_no_ct_value <- ipw.log.rank(times = data_no_ct$delay_os, failures = data_no_ct$status_vital,
									variable = data_no_ct$med_numeric, weights =data_no_ct$IPW_logistic_stab_trunc_right)
	print(pB_no_ct_value)
	
	survB_no_ctKM_weighted_df <- survB_no_ctKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
	head(survB_no_ctKM_weighted_df)	
		
	survB_no_ctKM_weighted_df$strata <- factor(survB_no_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
	p_val_toprint <- ifelse(pB_no_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pB_no_ct_value$p.value,3))))
	print(p_val_toprint)
	
	survB_no_ct_plot <- ggsurvplot(fit =survB_no_ctKM_weighted_df,
							palette = c(atc_color,"#868686"),
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
	survB_no_ct_plot <- survB_no_ct_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("...and by CT status"))),
							   atop(.(label_interaction))))) +
				labs(subtitle = label_no_ct) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold", color = "white")) 
	survB_no_ct_plot 
	
	#Save all plots in png
	plot_2G <- plot_grid(survA_plot + theme(plot.title = element_blank()), plotA$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
	save(survA,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survA.RData"))
	save(survAKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survAKM_weighted_df.RData"))
	png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page2G.png"),
		height = 2500, width = 1675, units = "px", res = 300,
		)
		print(plot_2G)
	dev.off() 
	
	plot_2H <- plot_grid(survB_ct_plot + theme(plot.title = element_blank()), plotB_ct$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
	save(survB_ct,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_ct.RData"))
	save(survB_ctKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_ctKM_weighted_df.RData"))
	png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page2H.png"),
		height = 2500, width = 1675, units = "px", res = 300,
		)
		print(plot_2H)
	dev.off() 
	
	plot_2I <- plot_grid(survB_no_ct_plot + theme(plot.title = element_blank()), plotB_no_ct$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
	save(survB_no_ct,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_no_ct.RData"))
	save(survB_no_ctKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_no_ctKM_weighted_df.RData"))
	png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page2I.png"),
		height = 2500, width = 1675, units = "px", res = 300,
		)
		print(plot_2I)
	dev.off() 
	
	#On assemble tous les plots
	row1 <- plot_grid(survA_plot, survB_ct_plot, survB_no_ct_plot,
						  plotA$table + theme(legend.position = "none"), 
						  plotB_ct$table + theme(legend.position = "none"), 
						  plotB_no_ct$table + theme(legend.position = "none"), 
						  ncol = 3 ,nrow = 2, rel_heights = c(3,1),
						  labels = c(label1,label2,'','','',''))
								                                       
	row1
						  
	return(row1)		

} 

################################################################################
#Kaplan-Meier survival curves endocrine therapy / no endocrine therapy
################################################################################

plot_survival_endocrine_therapy <- function(id,atc,atc_name,atc_color,
						 data_ct, data_no_ct, label1 = 'D', label2 = 'E',atc_level ="atc2",estimate  ="ATE",
						 ylab_txt = "Overall survival probability"){
		
	print(dim(id))
	text_all <- case_when(
		estimate == "ATT"  ~"Treated population",
		estimate == "matched" ~ "Matched population",
		TRUE ~ "Whole population"
	)
	label_facet = paste0(text_all, ", n = ",format(nrow(id), big.mark = " ")," (100%)")

	#Part 3 : plot B : survival - CT Yes/No

		#Test interaction CT vs no CT
		mod1 <- coxph(Surv(delay_os, status_vital) ~  ht + var + ht*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
		mod2 <- coxph(Surv(delay_os, status_vital) ~ ht + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
		tmp_anova <- anova(mod1,mod2,test = "Chisq") 
		label_interaction <- paste0("Interaction ",atc,"-HT p-value : ", 
								ifelse(round(tmp_anova[,"P(>|Chi|)"][2],3) == 0, "<0.001", round(tmp_anova[,"P(>|Chi|)"][2],3))) 
		print(label_interaction)
  
 			#plots
	 	print("Here just before row ct")
		label_ct <- paste0("Endocrine therapy, n = ", format(nrow(data_ct),big.mark = " "), " (", round(100*nrow(data_ct)/nrow(id),1),"%)")
		print(label_ct)
		print(head(data_ct))
		#Surv raw for risk table
		survB_ct <- survfit(Surv(delay_os, status_vital) ~ var, data = data_ct)
		pval_toprint_ct = surv_pvalue(survB_ct)$pval.txt
		plotB_ct <- ggsurvplot(fit = survB_ct, 
								data = data_ct,
								risk.table = TRUE, 
								pval = pval_toprint_ct, 
								ggtheme = theme_bw(), 
								palette = c(atc_color,"#868686"),
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = c("Yes","No"),
								censor = FALSE,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5) 
					
		#Survfit KM weighted
	survB_ctKM_weighted <- ipw.survival(times = data_ct$delay_os, failures = data_ct$status_vital,
									variable = data_ct$med_numeric, weights = data_ct$IPW_logistic_stab_trunc_right)
								
	pB_ct_value <- ipw.log.rank(times = data_ct$delay_os, failures = data_ct$status_vital,
									variable = data_ct$med_numeric, weights = data_ct$IPW_logistic_stab_trunc_right)
	print(pB_ct_value)
	
	survB_ctKM_weighted_df <- survB_ctKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
	head(survB_ctKM_weighted_df)	
		
	survB_ctKM_weighted_df$strata <- factor(survB_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
	p_val_toprint <- ifelse(pB_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pB_ct_value$p.value,3))))
	print(p_val_toprint)
	
	survB_ct_plot <- ggsurvplot(fit = survB_ctKM_weighted_df,
							palette = c(atc_color,"#868686"),
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
	survB_ct_plot <- survB_ct_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("Mortality by ",ifelse(nchar(atc_name) > 20, atc,atc_name), " and CT status"))),
							   atop(.(label_interaction))))) +
				labs(subtitle = label_ct) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold"))
	survB_ct_plot  
	
	#No CT
	label_no_ct <- paste0("No endocrine therapy, n = ", format(nrow(data_no_ct),big.mark = " "), " (", round(100*nrow(data_no_ct)/nrow(id),1),"%)")
	#Surv raw for risk table
		survB_no_ct <- survfit(Surv(delay_os, status_vital) ~ var, data = data_no_ct)
		pval_toprint_no_ct = surv_pvalue(survB_no_ct)$pval.txt
		plotB_no_ct <- ggsurvplot(fit = survB_no_ct, 
								data = data_ct,
								risk.table = TRUE, 
								pval = pval_toprint_no_ct, 
								ggtheme = theme_bw(), 
								palette = c(atc_color,"#868686"),
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = c("Yes","No"),
								censor = FALSE,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5) 
	
			#Survfit KM weighted
	survB_no_ctKM_weighted <- ipw.survival(times = data_no_ct$delay_os, failures = data_no_ct$status_vital,
									variable = data_no_ct$med_numeric, weights = data_no_ct$IPW_logistic_stab_trunc_right)
								
	pB_no_ct_value <- ipw.log.rank(times = data_no_ct$delay_os, failures = data_no_ct$status_vital,
									variable = data_no_ct$med_numeric, weights =data_no_ct$IPW_logistic_stab_trunc_right)
	print(pB_no_ct_value)
	
	survB_no_ctKM_weighted_df <- survB_no_ctKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
	head(survB_no_ctKM_weighted_df)	
		
	survB_no_ctKM_weighted_df$strata <- factor(survB_no_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
	p_val_toprint <- ifelse(pB_no_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pB_no_ct_value$p.value,3))))
	print(p_val_toprint)
	
	survB_no_ct_plot <- ggsurvplot(fit =survB_no_ctKM_weighted_df,
							palette = c(atc_color,"#868686"),
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
	survB_no_ct_plot <- survB_no_ct_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("...and by CT status"))),
							   atop(.(label_interaction))))) +
				labs(subtitle = label_no_ct) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold", color = "white")) 
	survB_no_ct_plot 

	#Save all plots in png	
	plot_2H <- plot_grid(survB_ct_plot + theme(plot.title = element_blank()), plotB_ct$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
	save(survB_ct,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_nplus.RData"))
	save(survB_ctKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_nplusKM_weighted_df.RData"))
	png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page2Hter.png"),
		height = 2500, width = 1675, units = "px", res = 300,
		)
		print(plot_2H)
	dev.off() 
	
	plot_2I <- plot_grid(survB_no_ct_plot + theme(plot.title = element_blank()), plotB_no_ct$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
	save(survB_no_ct,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_nmoins.RData"))
	save(survB_no_ctKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_nmoinsKM_weighted_df.RData"))
	png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page2Iter.png"),
		height = 2500, width = 1675, units = "px", res = 300,
		)
		print(plot_2I)
	dev.off() 
						  
	return(NULL)		

} 

################################################################################
#Kaplan-Meier survival curves node-positive / node-negative
################################################################################

plot_survival_node_status <- function(id,atc,atc_name,atc_color, data_ct, data_no_ct, label1 = 'D', label2 = 'E',atc_level ="atc2",estimate  ="ATE",
 ylab_txt = "Overall survival probability"){
		
	print(dim(id))
	text_all <- case_when(
		estimate == "ATT"  ~"Treated population",
		estimate == "matched" ~ "Matched population",
		TRUE ~ "Whole population"
	)
	label_facet = paste0(text_all, ", n = ",format(nrow(id), big.mark = " ")," (100%)")

#Part 3 : plot B : survival - CT Yes/No

		#Test interaction CT vs no CT
		mod1 <- coxph(Surv(delay_os, status_vital) ~  pnuicc_2cl + var + pnuicc_2cl*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
		mod2 <- coxph(Surv(delay_os, status_vital) ~ pnuicc_2cl + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
		tmp_anova <- anova(mod1,mod2,test = "Chisq") 
		label_interaction <- paste0("Interaction ",atc,"-CT p-value : ", 
								ifelse(round(tmp_anova[,"P(>|Chi|)"][2],3) == 0, "<0.001", round(tmp_anova[,"P(>|Chi|)"][2],3))) 
		print(label_interaction)
  
 			#plots
		label_ct <- paste0("Node-positive, n = ", format(nrow(data_ct),big.mark = " "), " (", round(100*nrow(data_ct)/nrow(id),1),"%)")
		#Surv raw for risk table
		survB_ct <- survfit(Surv(delay_os, status_vital) ~ var, data = data_ct)
		pval_toprint_ct = surv_pvalue(survB_ct)$pval.txt
		plotB_ct <- ggsurvplot(fit = survB_ct, 
								data = data_ct,
								risk.table = TRUE, 
								pval = pval_toprint_ct, 
								ggtheme = theme_bw(), 
								palette = c(atc_color,"#868686"),
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = c("Yes","No"),
								censor = FALSE,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5) 
					
		#Survfit KM weighted
	survB_ctKM_weighted <- ipw.survival(times = data_ct$delay_os, failures = data_ct$status_vital,
									variable = data_ct$med_numeric, weights = data_ct$IPW_logistic_stab_trunc_right)
								
	pB_ct_value <- ipw.log.rank(times = data_ct$delay_os, failures = data_ct$status_vital,
									variable = data_ct$med_numeric, weights = data_ct$IPW_logistic_stab_trunc_right)
	print(pB_ct_value)
	
	survB_ctKM_weighted_df <- survB_ctKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
	head(survB_ctKM_weighted_df)	
		
	survB_ctKM_weighted_df$strata <- factor(survB_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
	p_val_toprint <- ifelse(pB_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pB_ct_value$p.value,3))))
	print(p_val_toprint)
	
	survB_ct_plot <- ggsurvplot(fit = survB_ctKM_weighted_df,
							palette = c(atc_color,"#868686"),
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
	survB_ct_plot <- survB_ct_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("Mortality by ",ifelse(nchar(atc_name) > 20, atc,atc_name), " and CT status"))),
							   atop(.(label_interaction))))) +
				labs(subtitle = label_ct) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold"))
	survB_ct_plot  
	
	#No CT
	label_no_ct <- paste0("Node-negative, n = ", format(nrow(data_no_ct),big.mark = " "), " (", round(100*nrow(data_no_ct)/nrow(id),1),"%)")
	#Surv raw for risk table
		survB_no_ct <- survfit(Surv(delay_os, status_vital) ~ var, data = data_no_ct)
		pval_toprint_no_ct = surv_pvalue(survB_no_ct)$pval.txt
		plotB_no_ct <- ggsurvplot(fit = survB_no_ct, 
								data = data_ct,
								risk.table = TRUE, 
								pval = pval_toprint_no_ct, 
								ggtheme = theme_bw(), 
								palette = c(atc_color,"#868686"),
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = c("Yes","No"),
								censor = FALSE,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5) 
	
			#Survfit KM weighted
	survB_no_ctKM_weighted <- ipw.survival(times = data_no_ct$delay_os, failures = data_no_ct$status_vital,
									variable = data_no_ct$med_numeric, weights = data_no_ct$IPW_logistic_stab_trunc_right)
								
	pB_no_ct_value <- ipw.log.rank(times = data_no_ct$delay_os, failures = data_no_ct$status_vital,
									variable = data_no_ct$med_numeric, weights =data_no_ct$IPW_logistic_stab_trunc_right)
	print(pB_no_ct_value)
	
	survB_no_ctKM_weighted_df <- survB_no_ctKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
	head(survB_no_ctKM_weighted_df)	
		
	survB_no_ctKM_weighted_df$strata <- factor(survB_no_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
	p_val_toprint <- ifelse(pB_no_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pB_no_ct_value$p.value,3))))
	print(p_val_toprint)
	
	survB_no_ct_plot <- ggsurvplot(fit =survB_no_ctKM_weighted_df,
							palette = c(atc_color,"#868686"),
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
	survB_no_ct_plot <- survB_no_ct_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("...and by CT status"))),
							   atop(.(label_interaction))))) +
				labs(subtitle = label_no_ct) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold", color = "white")) 
	survB_no_ct_plot 

	#Save all plots in png	
	plot_2H <- plot_grid(survB_ct_plot + theme(plot.title = element_blank()), plotB_ct$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
	save(survB_ct,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_ht.RData"))
	save(survB_ctKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_htKM_weighted_df.RData"))
	png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page2Hbis.png"),
		height = 2500, width = 1675, units = "px", res = 300,
		)
		print(plot_2H)
	dev.off() 
	
	plot_2I <- plot_grid(survB_no_ct_plot + theme(plot.title = element_blank()), plotB_no_ct$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
	save(survB_no_ct,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_no_ht.RData"))
	save(survB_no_ctKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_no_htKM_weighted_df.RData"))
	png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page2Ibis.png"),
		height = 2500, width = 1675, units = "px", res = 300,
		)
		print(plot_2I)
	dev.off() 
						  
	return(NULL)		

}

################################################################################
#Kaplan-Meier survival curves by subtype
################################################################################

plot_survival_subtype <- function(id,id_subtype,atc,atc_name,atc_color,data_luminal, data_her2,data_tnbc, atc_level = "atc2",
estimate = "ATE",ylab_txt = "Overall survival probability"){

		#Test interaction subtype
		mod1 <- coxph(Surv(delay_os, status_vital) ~  subtype + var + subtype*var, data = id_subtype, weights = id_subtype$IPW_logistic_stab_trunc_right, robust = TRUE)
		mod2 <- coxph(Surv(delay_os, status_vital) ~  subtype + var, data = id_subtype, weights = id_subtype$IPW_logistic_stab_trunc_right, robust = TRUE)
		tmp_anova <- anova(mod1,mod2,test = "Chisq") 
		label_interaction <- paste0("Interaction ",
						atc,
						"-subtype p-value : ", 
						ifelse(round(tmp_anova[,"P(>|Chi|)"][2],3) == 0, "<0.001", round(tmp_anova[,"P(>|Chi|)"][2],3))
						) 
		label_interaction
		
		
		label_luminal <- paste0("Luminal, n = ", format(nrow(data_luminal),big.mark = " "), " (", round(100*nrow(data_luminal)/nrow(id_subtype),1),"%)")
		#Raw  data for risk table
		survC_luminal <- survfit(Surv(delay_os, status_vital) ~ var, data = data_luminal)
		pval_toprint_luminal = surv_pvalue(survC_luminal)$pval.txt
		legend_labs <- case_when(
			all(as.character(data_luminal$var) ==  "Yes") ~ list(c("Yes")),
			all(as.character(data_luminal$var) ==  "No") ~  list(c("No")),
			TRUE ~ list(c("Yes","No")),
		)
		palette_labs <- case_when(
			all(as.character(data_luminal$var) ==  "Yes") ~ list(c(atc_color)),
			all(as.character(data_luminal$var) ==  "No") ~  list(c("#868686")),
			TRUE ~ list(c(atc_color,"#868686")),
		)
		plotC_luminal <- ggsurvplot(fit = survC_luminal,
							data = id,
							risk.table = TRUE,
							pval = pval_toprint_luminal,
							ggtheme = theme_bw(), 
							palette = palette_labs[[1]],
							risk.table.col = "strata",
							risk.table.title = "Number at risk",
							legend = "none",
							legend.title = "",
							legend.labs = legend_labs[[1]],
							censor = F,
							ylab = ylab_txt,
							risk.table.fontsize = 3.5)
							
		#Weighted
		survC_luminalKM_weighted <- ipw.survival(times = data_luminal$delay_os, failures = data_luminal$status_vital,
							variable = data_luminal$med_numeric, weights = data_luminal$IPW_logistic_stab_trunc_right)
								
		pC_luminal_value <- ipw.log.rank(times = data_luminal$delay_os, failures = data_luminal$status_vital,
									variable = data_luminal$med_numeric, weights = data_luminal$IPW_logistic_stab_trunc_right)
		print(pC_luminal_value)
		
		survC_luminalKM_weighted_df <- survC_luminalKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
		head(survC_luminalKM_weighted_df)	
		
		survC_luminalKM_weighted_df$strata <- factor(survC_luminalKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
		p_val_toprint <- ifelse(pC_luminal_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pC_luminal_value$p.value,3))))
		print(p_val_toprint)
		
		survC_luminal_plot <- ggsurvplot(fit = survC_luminalKM_weighted_df,
							palette = palette_labs[[1]],
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
		survC_luminal_plot <- survC_luminal_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("Mortality by ",ifelse(nchar(atc_name) > 20, atc,atc_name), " and subtype"))),
								   atop(.(label_interaction))))) +
					labs(subtitle = label_luminal) +
					theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
					theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold"))
		survC_luminal_plot 
		
		#HER2+
		label_her2 <- paste0("HER2+, n = ", format(nrow(data_her2),big.mark = " "), " (", round(100*nrow(data_her2)/nrow(id_subtype),1),"%)")
		#Raw  data for risk table
		survC_her2 <- survfit(Surv(delay_os, status_vital) ~ var, data = data_her2)
		pval_toprint_her2 = surv_pvalue(survC_her2)$pval.txt
		legend_labs <- case_when(
			all(as.character(data_her2$var) ==  "Yes") ~ list(c("Yes")),
			all(as.character(data_her2$var) ==  "No") ~  list(c("No")),
			TRUE ~ list(c("Yes","No")),
		)
		palette_labs <- case_when(
			all(as.character(data_her2$var) ==  "Yes") ~ list(c(atc_color)),
			all(as.character(data_her2$var) ==  "No") ~  list(c("#868686")),
			TRUE ~ list(c(atc_color,"#868686")),
		)
		plotC_her2 <- ggsurvplot(fit = survC_her2,
							data = id,
							risk.table = TRUE,
							pval = pval_toprint_her2,
							ggtheme = theme_bw(), 
							palette = palette_labs[[1]],
							risk.table.col = "strata",
							risk.table.title = "Number at risk",
							legend = "none",
							legend.title = "",
							legend.labs = legend_labs[[1]],
							censor = F,
							ylab = ylab_txt,
							risk.table.fontsize = 3.5)
							
		#Weighted
		survC_her2KM_weighted <- ipw.survival(times = data_her2$delay_os, failures = data_her2$status_vital,
							variable = data_her2$med_numeric, weights = data_her2$IPW_logistic_stab_trunc_right)
								
		pC_her2_value <- ipw.log.rank(times = data_her2$delay_os, failures = data_her2$status_vital,
									variable = data_her2$med_numeric, weights = data_her2$IPW_logistic_stab_trunc_right)
		print(pC_her2_value)
		
		survC_her2KM_weighted_df <- survC_her2KM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
		head(survC_her2KM_weighted_df)	
		
		survC_her2KM_weighted_df$strata <- factor(survC_her2KM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
		p_val_toprint <- ifelse(pC_her2_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pC_her2_value$p.value,3))))
		print(p_val_toprint)
		
		survC_her2_plot <- ggsurvplot(fit = survC_her2KM_weighted_df,
							palette = palette_labs[[1]],
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
		survC_her2_plot <- survC_her2_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			labs(subtitle = label_her2) +
			theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))
		survC_her2_plot
		
		#TNBC
		label_tnbc <- paste0("TNBC, n = ", format(nrow(data_tnbc),big.mark = " "), " (", round(100*nrow(data_tnbc)/nrow(id_subtype),1),"%)")
		#Raw  data for risk table
		survC_tnbc <- survfit(Surv(delay_os, status_vital) ~ var, data = data_tnbc)
		pval_toprint_tnbc = surv_pvalue(survC_tnbc)$pval.txt
		legend_labs <- case_when(
			all(as.character(data_tnbc$var) ==  "Yes") ~ list(c("Yes")),
			all(as.character(data_tnbc$var) ==  "No") ~  list(c("No")),
			TRUE ~ list(c("Yes","No")),
		)
		palette_labs <- case_when(
			all(as.character(data_tnbc$var) ==  "Yes") ~ list(c(atc_color)),
			all(as.character(data_tnbc$var) ==  "No") ~  list(c("#868686")),
			TRUE ~ list(c(atc_color,"#868686")),
		)
		plotC_tnbc <- ggsurvplot(fit = survC_tnbc,
							data = id,
							risk.table = TRUE,
							pval = pval_toprint_tnbc,
							ggtheme = theme_bw(), 
							palette = palette_labs[[1]],
							risk.table.col = "strata",
							risk.table.title = "Number at risk",
							legend = "none",
							legend.title = "",
							legend.labs = legend_labs[[1]],
							censor = F,
							ylab = ylab_txt,
							risk.table.fontsize = 3.5)
							
		#Weighted
		survC_tnbcKM_weighted <- ipw.survival(times = data_tnbc$delay_os, failures = data_tnbc$status_vital,
							variable = data_tnbc$med_numeric, weights = data_tnbc$IPW_logistic_stab_trunc_right)
								
		pC_tnbc_value <- ipw.log.rank(times = data_tnbc$delay_os, failures = data_tnbc$status_vital,
									variable = data_tnbc$med_numeric, weights = data_tnbc$IPW_logistic_stab_trunc_right)
		print(pC_tnbc_value)
		
		survC_tnbcKM_weighted_df <- survC_tnbcKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
		head(survC_tnbcKM_weighted_df)	
		
		survC_tnbcKM_weighted_df$strata <- factor(survC_tnbcKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
		p_val_toprint <- ifelse(pC_tnbc_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pC_tnbc_value$p.value,3))))
		print(p_val_toprint)
		
		survC_tnbc_plot <- ggsurvplot(fit = survC_tnbcKM_weighted_df,
							palette = palette_labs[[1]],
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
		survC_tnbc_plot <- survC_tnbc_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			labs(subtitle = label_tnbc) +
			theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5))) 
		survC_tnbc_plot
		
			#Save all plots in png
		plot_3A <- plot_grid(survC_luminal_plot + theme(plot.title = element_blank()),plotC_luminal$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
		save(survC_luminal,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survC_luminal.RData"))
		save(survC_luminalKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survC_luminalKM_weighted_df.RData"))
		png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page3A.png"),
			height = 2500, width = 1675, units = "px", res = 300,
			)
			print(plot_3A)
		dev.off() 
		
		plot_3D <- plot_grid(survC_her2_plot  + theme(plot.title = element_blank()), plotC_her2$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
		save(survC_her2,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survC_her2.RData"))
		save(survC_her2KM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survC_her2KM_weighted_df.RData"))
		
		png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page3D.png"),
			height = 2500, width = 1675, units = "px", res = 300,
			)
			print(plot_3D)
		dev.off() 
		
		plot_3G <- plot_grid(survC_tnbc_plot + theme(plot.title = element_blank()), plotC_tnbc$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
		save(survC_tnbc,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survC_tnbc.RData"))
		save(survC_tnbcKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survC_tnbcKM_weighted_df.RData"))
		
		png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page3G.png"),
			height = 2500, width = 1675, units = "px", res = 300,
			)
			print(plot_3G)
		dev.off()  
		
		#On assemble tous les plots
		row234_col1 <- plot_grid(survC_luminal_plot, 
								plotC_luminal$table + theme(legend.position = "none"),
								survC_tnbc_plot, plotC_tnbc$table + theme(legend.position = "none"),
								survC_her2_plot, plotC_her2$table + theme(legend.position = "none"),							
						  nrow = 6 ,ncol = 1, rel_heights = c(3.35,1,3,1,3,1))

		 			
		return(row234_col1)
		
}  

################################################################################
# Kaplan-Meier survival curves by subtype in patients with chemotherapy
################################################################################

plot_survival_subtype_ct <- function(id,id_subtype,atc,atc_name,atc_color,data_ct_luminal, 
									   data_ct_her2,data_ct_tnbc, atc_level = "atc2", estimate = "ATE",
									    ylab_txt = "Overall survival probability"){

		#Test interaction CT vs no CT
		mod1 <- coxph(Surv(delay_os, status_vital) ~  subtype + var + ct + subtype*var*ct, data = id_subtype)
		mod2 <- coxph(Surv(delay_os, status_vital) ~ subtype + var + ct, data = id_subtype)
		tmp_anova <- anova(mod1,mod2,test = "Chisq") 
		label_interaction <- paste0("Interaction ",atc,"-subtype-CT p-value : ", 
									ifelse(round(tmp_anova[,"P(>|Chi|)"][2],3) == 0, "<0.001", round(tmp_anova[,"P(>|Chi|)"][2],3))) 
		label_interaction
		
		
		
		#Luminal - CT
		label_ct_luminal <- paste0("Chemotherapy- Luminal, n = ", format(nrow(data_ct_luminal),big.mark = " "), " (", round(100*nrow(data_ct_luminal)/nrow(id_subtype),1),"%)")
		survD_ct_luminal <- survfit(Surv(delay_os, status_vital) ~ var, data = data_ct_luminal)
		pval_toprint_ct_luminal = surv_pvalue(survD_ct_luminal)$pval.txt
		legend_labs <- case_when(
			all(as.character(data_ct_luminal$var) ==  "Yes") ~ list(c("Yes")),
			all(as.character(data_ct_luminal$var) ==  "No") ~  list(c("No")),
			TRUE ~ list(c("Yes","No")),
		)
		palette_labs <- case_when(
			all(as.character(data_ct_luminal$var) ==  "Yes") ~ list(c(atc_color)),
			all(as.character(data_ct_luminal$var) ==  "No") ~  list(c("#868686")),
			TRUE ~ list(c(atc_color,"#868686")),
		)
		#Raw
		plotD_ct_luminal <- ggsurvplot(fit = survD_ct_luminal,
							data = data_ct_luminal,
							risk.table = TRUE,
							pval = pval_toprint_ct_luminal,
							ggtheme = theme_bw(), 
							palette = palette_labs[[1]],
							risk.table.col = "strata",
							risk.table.title = "Number at risk",
							legend = "none",
							legend.title = "",
							legend.labs = legend_labs[[1]],
							censor = F,
							ylab = ylab_txt,
							risk.table.fontsize = 3.5)
							
		#Weighted KM
		survD_luminal_ctKM_weighted <- ipw.survival(times = data_ct_luminal$delay_os, failures = data_ct_luminal$status_vital,
							variable = data_ct_luminal$med_numeric, weights = data_ct_luminal$IPW_logistic_stab_trunc_right)
								
		pC_luminal_ct_value <- ipw.log.rank(times = data_ct_luminal$delay_os, failures = data_ct_luminal$status_vital,
									variable = data_ct_luminal$med_numeric, weights = data_ct_luminal$IPW_logistic_stab_trunc_right)
		print(pC_luminal_ct_value)
		
		survD_luminal_ctKM_weighted_df <- survD_luminal_ctKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
		head(survD_luminal_ctKM_weighted_df)	
		
		survD_luminal_ctKM_weighted_df$strata <- factor(survD_luminal_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
		p_val_toprint <- ifelse(pC_luminal_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pC_luminal_ct_value$p.value,3))))
		print(p_val_toprint)
		
		survD_luminal_ctplot <- ggsurvplot(fit = survD_luminal_ctKM_weighted_df,
							palette = palette_labs[[1]],
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
		survD_luminal_ctplot <- survD_luminal_ctplot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("Mortality by ",ifelse(nchar(atc_name) > 20, atc,atc_name), ", subtype and CT status"))),
						   atop(.(label_interaction))))) +
			labs(subtitle = label_ct_luminal) +
			theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
			theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold"))		
		
		#print(survD_luminal_ctplot)
		
		
		#her2 - CT
		label_ct_her2 <- paste0("Chemotherapy- HER2+, n = ", format(nrow(data_ct_her2),big.mark = " "), " (", round(100*nrow(data_ct_her2)/nrow(id_subtype),1),"%)")
		survD_ct_her2 <- survfit(Surv(delay_os, status_vital) ~ var, data = data_ct_her2)
		pval_toprint_ct_her2 = surv_pvalue(survD_ct_her2)$pval.txt
		legend_labs <- case_when(
			all(as.character(data_ct_her2$var) ==  "Yes") ~ list(c("Yes")),
			all(as.character(data_ct_her2$var) ==  "No") ~  list(c("No")),
			TRUE ~ list(c("Yes","No")),
		)
		palette_labs <- case_when(
			all(as.character(data_ct_her2$var) ==  "Yes") ~ list(c(atc_color)),
			all(as.character(data_ct_her2$var) ==  "No") ~  list(c("#868686")),
			TRUE ~ list(c(atc_color,"#868686")),
		)
		#Raw
		plotD_ct_her2 <- ggsurvplot(fit = survD_ct_her2,
							data = data_ct_her2,
							risk.table = TRUE,
							pval = pval_toprint_ct_her2,
							ggtheme = theme_bw(), 
							palette = palette_labs[[1]],
							risk.table.col = "strata",
							risk.table.title = "Number at risk",
							legend = "none",
							legend.title = "",
							legend.labs = legend_labs[[1]],
							censor = F,
							ylab = ylab_txt,
							risk.table.fontsize = 3.5)
							
		#Weighted KM
		survD_her2_ctKM_weighted <- ipw.survival(times = data_ct_her2$delay_os, failures = data_ct_her2$status_vital,
							variable = data_ct_her2$med_numeric, weights = data_ct_her2$IPW_logistic_stab_trunc_right)
								
		pC_her2_ct_value <- ipw.log.rank(times = data_ct_her2$delay_os, failures = data_ct_her2$status_vital,
									variable = data_ct_her2$med_numeric, weights = data_ct_her2$IPW_logistic_stab_trunc_right)
		print(pC_her2_ct_value)
		
		survD_her2_ctKM_weighted_df <- survD_her2_ctKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
		head(survD_her2_ctKM_weighted_df)	
		
		survD_her2_ctKM_weighted_df$strata <- factor(survD_her2_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
		p_val_toprint <- ifelse(pC_her2_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pC_her2_ct_value$p.value,3))))
		print(p_val_toprint)
		
		survD_her2_ctplot <- ggsurvplot(fit = survD_her2_ctKM_weighted_df,
							palette = palette_labs[[1]],
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
		survD_her2_ctplot <- survD_her2_ctplot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			labs(subtitle = label_ct_her2) +
			theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5))) 	
		
		#print(survD_her2_ctplot)
		
		#tnbc - CT
		label_ct_tnbc <- paste0("Chemotherapy- TNBC, n = ", format(nrow(data_ct_tnbc),big.mark = " "), " (", round(100*nrow(data_ct_tnbc)/nrow(id_subtype),1),"%)")
		survD_ct_tnbc <- survfit(Surv(delay_os, status_vital) ~ var, data = data_ct_tnbc)
		pval_toprint_ct_tnbc = surv_pvalue(survD_ct_tnbc)$pval.txt
		legend_labs <- case_when(
			all(as.character(data_ct_tnbc$var) ==  "Yes") ~ list(c("Yes")),
			all(as.character(data_ct_tnbc$var) ==  "No") ~  list(c("No")),
			TRUE ~ list(c("Yes","No")),
		)
		palette_labs <- case_when(
			all(as.character(data_ct_tnbc$var) ==  "Yes") ~ list(c(atc_color)),
			all(as.character(data_ct_tnbc$var) ==  "No") ~  list(c("#868686")),
			TRUE ~ list(c(atc_color,"#868686")),
		)
		#Raw
		plotD_ct_tnbc <- ggsurvplot(fit = survD_ct_tnbc,
							data = data_ct_tnbc,
							risk.table = TRUE,
							pval = pval_toprint_ct_tnbc,
							ggtheme = theme_bw(), 
							palette = palette_labs[[1]],
							risk.table.col = "strata",
							risk.table.title = "Number at risk",
							legend = "none",
							legend.title = "",
							legend.labs = legend_labs[[1]],
							censor = F,
							ylab = ylab_txt,
							risk.table.fontsize = 3.5)
							
		#Weighted KM
		survD_tnbc_ctKM_weighted <- ipw.survival(times = data_ct_tnbc$delay_os, failures = data_ct_tnbc$status_vital,
							variable = data_ct_tnbc$med_numeric, weights = data_ct_tnbc$IPW_logistic_stab_trunc_right)
								
		pC_tnbc_ct_value <- ipw.log.rank(times = data_ct_tnbc$delay_os, failures = data_ct_tnbc$status_vital,
									variable = data_ct_tnbc$med_numeric, weights = data_ct_tnbc$IPW_logistic_stab_trunc_right)
		print(pC_tnbc_ct_value)
		
		survD_tnbc_ctKM_weighted_df <- survD_tnbc_ctKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
		head(survD_tnbc_ctKM_weighted_df)	
		
		survD_tnbc_ctKM_weighted_df$strata <- factor(survD_tnbc_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
		p_val_toprint <- ifelse(pC_tnbc_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pC_tnbc_ct_value$p.value,3))))
		print(p_val_toprint)
		
		survD_tnbc_ctplot <- ggsurvplot(fit = survD_tnbc_ctKM_weighted_df,
							palette = palette_labs[[1]],
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
		survD_tnbc_ctplot <- survD_tnbc_ctplot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			labs(subtitle = label_ct_tnbc) +
			theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5))) 
		
		#Save all plots in png
		plot_3B <- plot_grid(survD_luminal_ctplot + theme(plot.title = element_blank()), plotD_ct_luminal$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
		save(survD_ct_luminal,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_ct_luminal.RData"))
		save(survD_luminal_ctKM_weighted,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_luminal_ctKM_weighted.RData"))
		png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page3B.png"),
			height = 2500, width = 1675, units = "px", res = 300,
			)
			print(plot_3B)
		dev.off() 
		
		plot_3E <- plot_grid(survD_her2_ctplot + theme(plot.title = element_blank()), plotD_ct_her2$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
		save(survD_ct_her2,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_ct_her2.RData"))
		save(survD_her2_ctKM_weighted,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_her2_ctKM_weighted.RData"))
		
		png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page3E.png"),
			height = 2500, width = 1675, units = "px", res = 300,
			)
			print(plot_3E)
		dev.off() 
		
		plot_3H <- plot_grid(survD_tnbc_ctplot + theme(plot.title = element_blank()), plotD_ct_tnbc$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
		save(survD_ct_tnbc,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_ct_tnbc.RData"))
		save(survD_tnbc_ctKM_weighted,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_tnbc_ctKM_weighted.RData"))
		
		png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page3H.png"),
			height = 2500, width = 1675, units = "px", res = 300,
			)
			print(plot_3H)
		dev.off() 
				
		row234_col2 <- plot_grid(survD_luminal_ctplot, plotD_ct_luminal$table + theme(legend.position = "none"),
								survD_tnbc_ctplot, plotD_ct_tnbc$table + theme(legend.position = "none"),
								survD_her2_ctplot, plotD_ct_her2$table + theme(legend.position = "none"),
								
						  nrow = 6 ,ncol = 1, rel_heights = c(3.35,1,3,1,3,1))
				 		 			
		return(row234_col2)
		
} 

################################################################################
# Kaplan-Meier survival curves by subtype in patients without chemotherapy
################################################################################

plot_survival_subtype_no_ct <- function(id,id_subtype,atc,atc_name,atc_color,data_no_ct_luminal, 
									   data_no_ct_her2,data_no_ct_tnbc, atc_level = "atc2", estimate = "ATE",
									   ylab_txt = "Overall survival probability"){
		
			#Test interaction CT vs no CT
			mod1 <- coxph(Surv(delay_os, status_vital) ~  subtype + var + ct + subtype*var*ct, data = id_subtype)
			mod2 <- coxph(Surv(delay_os, status_vital) ~ subtype + var + ct, data = id_subtype)
			tmp_anova <- anova(mod1,mod2,test = "Chisq") 
			label_interaction <- paste0("Interaction ",atc,"-subtype-CT p-value : ", 
								ifelse(round(tmp_anova[,"P(>|Chi|)"][2],3) == 0, "<0.001", round(tmp_anova[,"P(>|Chi|)"][2],3))) 
			label_interaction
			
			#Luminal - no_ct
			label_no_ct_luminal <- paste0("No chemotherapy- Luminal, n = ", format(nrow(data_no_ct_luminal),big.mark = " "), " (", round(100*nrow(data_no_ct_luminal)/nrow(id_subtype),1),"%)")
			survD_no_ct_luminal <- survfit(Surv(delay_os, status_vital) ~ var, data = data_no_ct_luminal)
			pval_toprint_no_ct_luminal = surv_pvalue(survD_no_ct_luminal)$pval.txt
			legend_labs <- case_when(
				all(as.character(data_no_ct_luminal$var) ==  "Yes") ~ list(c("Yes")),
				all(as.character(data_no_ct_luminal$var) ==  "No") ~  list(c("No")),
				TRUE ~ list(c("Yes","No")),
			)
			palette_labs <- case_when(
				all(as.character(data_no_ct_luminal$var) ==  "Yes") ~ list(c(atc_color)),
				all(as.character(data_no_ct_luminal$var) ==  "No") ~  list(c("#868686")),
				TRUE ~ list(c(atc_color,"#868686")),
			)
			#Raw
			plotD_no_ct_luminal <- ggsurvplot(fit = survD_no_ct_luminal,
								data = data_no_ct_luminal,
								risk.table = TRUE,
								pval = pval_toprint_no_ct_luminal,
								ggtheme = theme_bw(), 
								palette = palette_labs[[1]],
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = legend_labs[[1]],
								censor = F,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5)
								
			#Weighted KM
			survD_luminal_no_ctKM_weighted <- ipw.survival(times = data_no_ct_luminal$delay_os, failures = data_no_ct_luminal$status_vital,
								variable = data_no_ct_luminal$med_numeric, weights = data_no_ct_luminal$IPW_logistic_stab_trunc_right)
									
			pC_luminal_no_ct_value <- ipw.log.rank(times = data_no_ct_luminal$delay_os, failures = data_no_ct_luminal$status_vital,
										variable = data_no_ct_luminal$med_numeric, weights = data_no_ct_luminal$IPW_logistic_stab_trunc_right)
			print(pC_luminal_no_ct_value)
			
			survD_luminal_no_ctKM_weighted_df <-survD_luminal_no_ctKM_weighted$table.surv %>%
			as.data.frame() %>%
			rename(time = times, surv = survival, strata = variable) %>%
			mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
			mutate(upper = surv) %>%
			mutate(lower = surv)
			head(survD_luminal_no_ctKM_weighted_df)	
			
			survD_luminal_no_ctKM_weighted_df$strata <- factor(survD_luminal_no_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
			p_val_toprint <- ifelse(pC_luminal_no_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pC_luminal_no_ct_value$p.value,3))))
			print(p_val_toprint)
			
			survD_luminal_no_ctplot <- ggsurvplot(fit = survD_luminal_no_ctKM_weighted_df,
								palette = palette_labs[[1]],
								#pval = 0.03#,
								censor  = FALSE,
								legend.title = name_atc_tokeep,
								ggtheme = theme_bw(),
								ylab = ylab_txt,
								size = 1
								)
								
			survD_luminal_no_ctplot <- survD_luminal_no_ctplot + 
				annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
				ggtitle(bquote(atop(bold(.(paste0("Mortality by ",ifelse(nchar(atc_name) > 20, atc,atc_name), ", subtype and CT status"))),
							   atop(.(label_interaction))))) +
				labs(subtitle = label_no_ct_luminal) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold", color = "white"))		
			
			#print(survD_luminal_no_ctplot)
			
			label_no_ct_her2 <- paste0("No chemotherapy- HER2+, n = ", format(nrow(data_no_ct_her2),big.mark = " "), " (", round(100*nrow(data_no_ct_her2)/nrow(id_subtype),1),"%)")
			survD_no_ct_her2 <- survfit(Surv(delay_os, status_vital) ~ var, data = data_no_ct_her2)
			pval_toprint_no_ct_her2 = surv_pvalue(survD_no_ct_her2)$pval.txt
			legend_labs <- case_when(
				all(as.character(data_no_ct_her2$var) ==  "Yes") ~ list(c("Yes")),
				all(as.character(data_no_ct_her2$var) ==  "No") ~  list(c("No")),
				TRUE ~ list(c("Yes","No")),
			)
			palette_labs <- case_when(
				all(as.character(data_no_ct_her2$var) ==  "Yes") ~ list(c(atc_color)),
				all(as.character(data_no_ct_her2$var) ==  "No") ~  list(c("#868686")),
				TRUE ~ list(c(atc_color,"#868686")),
			)
			#Raw
			plotD_no_ct_her2 <- ggsurvplot(fit = survD_no_ct_her2,
								data = data_no_ct_her2,
								risk.table = TRUE,
								pval = pval_toprint_no_ct_her2,
								ggtheme = theme_bw(), 
								palette = palette_labs[[1]],
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = legend_labs[[1]],
								censor = F,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5)
								
			#Weighted KM
			survD_her2_no_ctKM_weighted <- ipw.survival(times = data_no_ct_her2$delay_os, failures = data_no_ct_her2$status_vital,
								variable = data_no_ct_her2$med_numeric, weights = data_no_ct_her2$IPW_logistic_stab_trunc_right)
									
			pC_her2_no_ct_value <- ipw.log.rank(times = data_no_ct_her2$delay_os, failures = data_no_ct_her2$status_vital,
										variable = data_no_ct_her2$med_numeric, weights = data_no_ct_her2$IPW_logistic_stab_trunc_right)
			print(pC_her2_no_ct_value)
			
			survD_her2_no_ctKM_weighted_df <- survD_her2_no_ctKM_weighted$table.surv %>%
			as.data.frame() %>%
			rename(time = times, surv = survival, strata = variable) %>%
			mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
			mutate(upper = surv) %>%
			mutate(lower = surv)
			head(survD_her2_no_ctKM_weighted_df)	
			
			survD_her2_no_ctKM_weighted_df$strata <- factor(survD_her2_no_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
			p_val_toprint <- ifelse(pC_her2_no_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pC_her2_no_ct_value$p.value,3))))
			print(p_val_toprint)
			
			survD_her2_no_ctplot <- ggsurvplot(fit = survD_her2_no_ctKM_weighted_df,
								palette = palette_labs[[1]],
								#pval = 0.03#,
								censor  = FALSE,
								legend.title = name_atc_tokeep,
								ggtheme = theme_bw(),
								ylab = ylab_txt,
								size = 1
								)
								
			survD_her2_no_ctplot <- survD_her2_no_ctplot + 
				annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
				labs(subtitle = label_no_ct_her2) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))	
			
			#print(survD_her2_no_ctplot)

			label_no_ct_undefined <- paste0("No chemotherapy- Undefined, n = ", format(nrow(data_no_ct_undefined),big.mark = " "), " (", round(100*nrow(data_no_ct_undefined)/nrow(id_subtype),1),"%)")
			survD_no_ct_undefined <- survfit(Surv(delay_os, status_vital) ~ var, data = data_no_ct_undefined)
			pval_toprint_no_ct_undefined = surv_pvalue(survD_no_ct_undefined)$pval.txt
			legend_labs <- case_when(
				all(as.character(data_no_ct_undefined$var) ==  "Yes") ~ list(c("Yes")),
				all(as.character(data_no_ct_undefined$var) ==  "No") ~  list(c("No")),
				TRUE ~ list(c("Yes","No")),
			)
			palette_labs <- case_when(
				all(as.character(data_no_ct_undefined$var) ==  "Yes") ~ list(c(atc_color)),
				all(as.character(data_no_ct_undefined$var) ==  "No") ~  list(c("#868686")),
				TRUE ~ list(c(atc_color,"#868686")),
			)
			#Raw
			plotD_no_ct_undefined <- ggsurvplot(fit = survD_no_ct_undefined,
								data = data_no_ct_undefined,
								risk.table = TRUE,
								pval = pval_toprint_no_ct_undefined,
								ggtheme = theme_bw(), 
								palette = palette_labs[[1]],
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = legend_labs[[1]],
								censor = F,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5)
								
			#Weighted KM
			survD_undefined_no_ctKM_weighted <- ipw.survival(times = data_no_ct_undefined$delay_os, failures = data_no_ct_undefined$status_vital,
								variable = data_no_ct_undefined$med_numeric, weights = data_no_ct_undefined$IPW_logistic_stab_trunc_right)
									
			pC_undefined_no_ct_value <- ipw.log.rank(times = data_no_ct_undefined$delay_os, failures = data_no_ct_undefined$status_vital,
										variable = data_no_ct_undefined$med_numeric, weights = data_no_ct_undefined$IPW_logistic_stab_trunc_right)
			print(pC_undefined_no_ct_value)
			
			survD_undefined_no_ctKM_weighted_df <- survD_undefined_no_ctKM_weighted$table.surv %>%
			as.data.frame() %>%
			rename(time = times, surv = survival, strata = variable) %>%
			mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
			mutate(upper = surv) %>%
			mutate(lower = surv)
			head(survD_undefined_no_ctKM_weighted_df)	
			
			survD_undefined_no_ctKM_weighted_df$strata <- factor(survD_undefined_no_ctKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
			p_val_toprint <- ifelse(pC_undefined_no_ct_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pC_undefined_no_ct_value$p.value,3))))
			print(p_val_toprint)
			
			survD_undefined_no_ctplot <- ggsurvplot(fit = survD_undefined_no_ctKM_weighted_df,
								palette = palette_labs[[1]],
								#pval = 0.03#,
								censor  = FALSE,
								legend.title = name_atc_tokeep,
								ggtheme = theme_bw(),
								ylab = ylab_txt,
								size = 1
								)
								
			survD_undefined_no_ctplot <- survD_undefined_no_ctplot + 
				annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
				labs(subtitle = label_no_ct_undefined) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  
			
			#print(survD_undefined_no_ctplot)
					#Save all plots in png
		plot_3C <- plot_grid(survD_luminal_no_ctplot + theme(plot.title = element_blank()), plotD_no_ct_luminal$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
		save(survD_no_ct_luminal,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_no_ct_luminal.RData"))
		save(survD_luminal_no_ctKM_weighted,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_luminal_no_ctKM_weighted.RData"))
		
		png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page3C.png"),
			height = 2500, width = 1675, units = "px", res = 300,
			)
			print(plot_3C)
		dev.off() 
		
		plot_3F <- plot_grid(survD_her2_no_ctplot + theme(plot.title = element_blank()), plotD_no_ct_her2$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
		save(survD_no_ct_her2,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_no_ct_her2.RData"))
		save(survD_her2_no_ctKM_weighted,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_her2_no_ctKM_weighted.RData"))
		
		png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page3F.png"),
			height = 2500, width = 1675, units = "px", res = 300,
			)
			print(plot_3F)
		dev.off() 
		
		plot_3I <- plot_grid(survD_undefined_no_ctplot + theme(plot.title = element_blank()), plotD_no_ct_undefined$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
		save(survD_no_ct_undefined,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_no_ct_undefined.RData"))
		save(survD_undefined_no_ctKM_weighted,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"survD_undefined_no_ctKM_weighted.RData"))
		
		png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page3I.png"),
			height = 2500, width = 1675, units = "px", res = 300,
			)
			print(plot_3I)
		dev.off() 
								
			#On assemble tous les plots
			row234_col3 <- plot_grid(survD_luminal_no_ctplot, plotD_no_ct_luminal$table + theme(legend.position = "none"),
									survD_undefined_no_ctplot, plotD_no_ct_undefined$table + theme(legend.position = "none"),
									survD_her2_no_ctplot, plotD_no_ct_her2$table + theme(legend.position = "none"),								
							  nrow = 6 ,ncol = 1, rel_heights = c(3.35,1,3,1,3,1))
							  
		 			
		return(row234_col3)
		
} 

################################################################################
# Kaplan-Meier survival curves by age classes
################################################################################

plot_survival_age <- function(id,atc,atc_name,atc_color,
									data_young,
									data_middle,
									data_old,
									label1 = 'D', label2 = 'E', label3 = 'F', 
									 atc_level ="atc2",estimate  ="ATE",
									  ylab_txt = "Overall survival probability"){
		
	print(dim(id))
	text_all <- case_when(
		estimate == "ATT"  ~"Treated population",
		estimate == "matched" ~ "Matched population",
		TRUE ~ "Whole population"
	)
	label_facet = paste0(text_all, ", n = ",format(nrow(id), big.mark = " ")," (100%)")

#Survival by age

		#Test interaction comed with age classes
		mod1 <- coxph(Surv(delay_os, status_vital) ~  age_3cl + var + age_3cl*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
		mod2 <- coxph(Surv(delay_os, status_vital) ~ age_3cl + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
		tmp_anova <- anova(mod1,mod2,test = "Chisq") 
		label_interaction <- paste0("Interaction ",atc,"-CT p-value : ", 
								ifelse(round(tmp_anova[,"P(>|Chi|)"][2],3) == 0, "<0.001", round(tmp_anova[,"P(>|Chi|)"][2],3))) 
		print(label_interaction)
  
 		#Plots
		label_age_young <- paste0("< 50 years old, n = ", format(nrow(data_young),big.mark = " "), " (", round(100*nrow(data_young)/nrow(id),1),"%)")
		#Surv raw for risk table
		survB_age_young <- survfit(Surv(delay_os, status_vital) ~ var, data = data_young)
		pval_toprint_age_young = surv_pvalue(survB_age_young)$pval.txt
		plotB_age_young <- ggsurvplot(fit = survB_age_young, 
								data = data_young,
								risk.table = TRUE, 
								pval = pval_toprint_age_young, 
								ggtheme = theme_bw(), 
								palette = c(atc_color,"#868686"),
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = c("Yes","No"),
								censor = FALSE,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5) 
					
		#Survfit KM weighted
	survB_age_youngKM_weighted <- ipw.survival(times = data_young$delay_os, failures = data_young$status_vital,
									variable = data_young$med_numeric, weights = data_young$IPW_logistic_stab_trunc_right)
								
	pB_age_young_value <- ipw.log.rank(times = data_young$delay_os, failures = data_young$status_vital,
									variable = data_young$med_numeric, weights = data_young$IPW_logistic_stab_trunc_right)
	print(pB_age_young_value)
	
	survB_age_youngKM_weighted_df <- survB_age_youngKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
	head(survB_age_youngKM_weighted_df)	
		
	survB_age_youngKM_weighted_df$strata <- factor(survB_age_youngKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
	p_val_toprint <- ifelse(pB_age_young_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pB_age_young_value$p.value,3))))
	print(p_val_toprint)
	
	survB_age_young_plot <- ggsurvplot(fit = survB_age_youngKM_weighted_df,
							palette = c(atc_color,"#868686"),
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
	survB_age_young_plot <- survB_age_young_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("Mortality by ",ifelse(nchar(atc_name) > 20, atc,atc_name), " and age classes"))),
							   atop(.(label_interaction))))) +
				labs(subtitle = label_age_young) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold"))
	survB_age_young_plot
	
	#Age middle
	labelage_middle <- paste0("50-75 years old, n = ", format(nrow(data_middle),big.mark = " "), " (", round(100*nrow(data_middle)/nrow(id),1),"%)")
	#Surv raw for risk table
		survB_age_middle <- survfit(Surv(delay_os, status_vital) ~ var, data = data_middle)
		pval_toprintage_middle = surv_pvalue(survB_age_middle)$pval.txt
		plotB_age_middle <- ggsurvplot(fit = survB_age_middle, 
								data = data_middle,
								risk.table = TRUE, 
								pval = pval_toprintage_middle, 
								ggtheme = theme_bw(), 
								palette = c(atc_color,"#868686"),
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = c("Yes","No"),
								censor = FALSE,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5) 
	
			#Survfit KM weighted
	survB_age_middleKM_weighted <- ipw.survival(times = data_middle$delay_os, failures = data_middle$status_vital,
									variable = data_middle$med_numeric, weights = data_middle$IPW_logistic_stab_trunc_right)
								
	pBage_middle_value <- ipw.log.rank(times = data_middle$delay_os, failures = data_middle$status_vital,
									variable = data_middle$med_numeric, weights =data_middle$IPW_logistic_stab_trunc_right)
	print(pBage_middle_value)
	
	survB_age_middleKM_weighted_df <- survB_age_middleKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
	head(survB_age_middleKM_weighted_df)	
		
	survB_age_middleKM_weighted_df$strata <- factor(survB_age_middleKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
	p_val_toprint <- ifelse(pBage_middle_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pBage_middle_value$p.value,3))))
	print(p_val_toprint)
	
	survB_age_middle_plot <- ggsurvplot(fit =survB_age_middleKM_weighted_df,
							palette = c(atc_color,"#868686"),
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
	survB_age_middle_plot <- survB_age_middle_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("...and by CT status"))),
							   atop(.(label_interaction))))) +
				labs(subtitle = labelage_middle) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold", color = "white")) 
	survB_age_middle_plot 
	
		#Age middle
	labelage_old <- paste0(">75 years old, n = ", format(nrow(data_old),big.mark = " "), " (", round(100*nrow(data_old)/nrow(id),1),"%)")
	#Surv raw for risk table
		survB_age_old <- survfit(Surv(delay_os, status_vital) ~ var, data = data_old)
		pval_toprintage_old = surv_pvalue(survB_age_old)$pval.txt
		plotB_age_old <- ggsurvplot(fit = survB_age_old, 
								data = data_old,
								risk.table = TRUE, 
								pval = pval_toprintage_old, 
								ggtheme = theme_bw(), 
								palette = c(atc_color,"#868686"),
								risk.table.col = "strata",
								risk.table.title = "Number at risk",
								legend = "none",
								legend.title = "",
								legend.labs = c("Yes","No"),
								censor = FALSE,
								ylab = ylab_txt,
								risk.table.fontsize = 3.5) 
	
			#Survfit KM weighted
	survB_age_oldKM_weighted <- ipw.survival(times = data_old$delay_os, failures = data_old$status_vital,
									variable = data_old$med_numeric, weights = data_old$IPW_logistic_stab_trunc_right)
								
	pBage_old_value <- ipw.log.rank(times = data_old$delay_os, failures = data_old$status_vital,
									variable = data_old$med_numeric, weights =data_old$IPW_logistic_stab_trunc_right)
	print(pBage_old_value)
	
	survB_age_oldKM_weighted_df <- survB_age_oldKM_weighted$table.surv %>%
		as.data.frame() %>%
		rename(time = times, surv = survival, strata = variable) %>%
		mutate(n.censor = n.risk - lead(n.risk) - n.event) %>%
		mutate(upper = surv) %>%
		mutate(lower = surv)
	head(survB_age_oldKM_weighted_df)	
		
	survB_age_oldKM_weighted_df$strata <- factor(survB_age_oldKM_weighted_df$strata, levels = c(1,0), labels = c("Yes","No"))
	p_val_toprint <- ifelse(pBage_old_value$p.value < 0.001, "p < 0.001", paste0("p = ",as.character(round(pBage_old_value$p.value,3))))
	print(p_val_toprint)
	
	survB_age_old_plot <- ggsurvplot(fit =survB_age_oldKM_weighted_df,
							palette = c(atc_color,"#868686"),
							#pval = 0.03#,
							censor  = FALSE,
							legend.title = name_atc_tokeep,
							ggtheme = theme_bw(),
							ylab = ylab_txt,
							size = 1
							)
							
	survB_age_old_plot <- survB_age_old_plot + 
			annotate("text", x = 10, y = 0.25, label = p_val_toprint, cex = 5, col = "black") +
			ggtitle(bquote(atop(bold(.(paste0("...and by CT status"))),
							   atop(.(label_interaction))))) +
				labs(subtitle = labelage_old) +
				theme(plot.subtitle = element_textbox(hjust = 0.5, margin = margin(t=5, b =5)))  + 
				theme(plot.title = element_text(size = rel(1.2), hjust = 0, face = "bold", color = "white")) 
	survB_age_old_plot
	
	

	#Save all plots in png	
	plot_2H <- plot_grid(survB_age_young_plot + theme(plot.title = element_blank()), plotB_age_young$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
	save(survB_age_young,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_ht.RData"))
	save(survB_age_youngKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_htKM_weighted_df.RData"))
	png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page2Hage.png"),
		height = 2500, width = 1675, units = "px", res = 300,
		)
		print(plot_2H)
	dev.off()
	
	plot_2I <- plot_grid(survB_age_middle_plot + theme(plot.title = element_blank()), plotB_age_middle$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
	save(survB_age_middle,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_no_ht.RData"))
	save(survB_age_middleKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_no_htKM_weighted_df.RData"))
	png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page2Iage.png"),
		height = 2500, width = 1675, units = "px", res = 300,
		)
		print(plot_2I)
	dev.off() 
	
	plot_2K <- plot_grid(survB_age_old_plot + theme(plot.title = element_blank()), plotB_age_old$table + theme(legend.position = "none"), nrow = 2, rel_heights = c(3,1))
	save(survB_age_old,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_no_ht.RData"))
	save(survB_age_oldKM_weighted_df,file = paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_survB_no_htKM_weighted_df.RData"))
	png(paste0(folder_txt,atc_level,"/",estimate,"/survival/",atc,"_page2Jage.png"),
		height = 2500, width = 1675, units = "px", res = 300,
		)
		print(plot_2K)
	dev.off() 
						  
	return(NULL)

}