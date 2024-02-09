/* ******************************************************************** */
/* Estimation of ATE
- Cox Hazard Ratio
- Adjusted Kaplan-Meier survival curves

For whole population and by subgroups :
	-By BC subtype
	-By endocrine therapy status
	-By chemotherapy status
	-By age
	-By nodal status
*/
/* ******************************************************************** */

/* ******************************************************************** */
/* Step  1 : descriptive plots */
/* ******************************************************************** */
OPTIONS SET=R_HOME "C:\Program Files\R\R-3.6.3";
proc R;
export data = fresh.fresh_preprocessed_red R=df; /* Database with all variables */
export data = fresh.top_atc5_snds R = list_comed; /* List of medication to test */
submit;

#Libraries
library(tidyverse)
source("v10_util.R")
list_comed <- list_comed
df_tot <- df

#Loop on all medication
for (i in 1:nrow(list_comed)){

	#ATC code and label of the molecules
	atc_tokeep <- as.character(list_comed$PHA_ATC_C07[i]) <- atc
	name_atc_tokeep = list_comed$PHA_ATC_LIB[i] <- atc_name
	print(atc_tokeep)
	print(name_atc_tokeep)

	# ******************************************************************** #
	# Database preprocessing for the specific medication #
	# ******************************************************************** #

	id <- prepare_df_atc5(df_tot, atc_tokeep)

	 # ******************************************************************** #
	 # Get the color for the medication (depends on ATC level 1 code)  #
	 # ******************************************************************** #

	color_palette <- c("#0DB14B","#CC004C","#CC004C","#F37021",
							  "#FCB711","#0089D0","#6460AA","#006D77")
	names(color_palette) <- c("A","B","C","H","M","N","R","Others")
	atc_color <- as.character(color_palette[ifelse(substr(atc_tokeep,1,1) %in% names(color_palette),
																	 substr(atc_tokeep,1,1),
																	 "Others")])
	print(atc_color)

	# ******************************************************************** #
	# Create plot legend #
	# ******************************************************************** #

	my_hist <- ggplot(data=id, aes(var,fill= var)) +
		geom_bar() +
		scale_fill_manual(values = c(atc_color,"#868686")) +
		labs(fill = name_atc_tokeep) +
		guides(fill = guide_legend(nrow = 1, title.position = "left"))
	legend <- cowplot::get_legend(my_hist)

	  # ******************************************************************** #
      # First row of plots :
      # - bar plot for comedication intake yes versus no
      # - violin plot for age by comedication intake
      #  -violin plot for richness index by comedication intake
      # -percentage of patients with at least one comorbid condition by comedication intake
      # -percentage of obsese patients with  condition by comedication intake

      # Plot are saved individually in the folder 'descriptive'
      # And merged all together in the folder 'plot_merged'
      # ******************************************************************** #

	row_a1 <- plot_row1_boxplot(id,atc,atc_name,atc_color,atc_level = "atc5")

	# ******************************************************************** #
	# Second, third, fourth, and fifth rows of plot
	# Bar plots for all comorbid conditions by comedication intake
	# ******************************************************************** #

	row2345 <- plot_row2345_boxplot(id,atc_name,atc_color,atc_level = "atc5")

	# ******************************************************************** #
	# Merge all rows and save merged plot #
	# ******************************************************************** #

	plot_total1 <- plot_grid(row_a1, row2345,legend, ncol = 1 ,nrow = 3, rel_heights = c(1.5,4,0.1))
	pdf(file = paste0(folder_txt,atc,"_page1.pdf"), width=15.7, height=24, onefile = T)
						print(plot_total1)
	dev.off()

}

endsubmit;
run;

/* ******************************************************************** */
/* Step  2 : table one by comedication intake */
/* ******************************************************************** */

OPTIONS SET=R_HOME "C:\Program Files\R\R-3.6.3";
proc R;
export data = fresh.fresh_preprocessed_red R=df;  /* Database with all variables */
export data = fresh.top_atc5_snds R = list_comed;  /* List of medication to test */
export data = fresh.ref_atc10 R = ref;  /* All ATC classes labels */
submit;

#Libraries
library(tidyverse)
source("v10_util.R")
#Pretty colnames for variables
pretty_colnames <- readxl::read_xlsx(paste0(folder_txt,"ref/pretty_colnames.xlsx"))
list_comed <- list_comed
df_tot <- df

#Loop on all medications to test
for (i in 1:nrow(list_comed)){

	#ATC code and label of the molecules
	atc_tokeep <- as.character(list_comed$PHA_ATC_C07[i]) <- atc
	name_atc_tokeep = list_comed$PHA_ATC_LIB[i] <- atc_name
	print(atc_tokeep)
	print(name_atc_tokeep)

	# ******************************************************************** #
	# Database preprocessing for the specific medication #
	# ******************************************************************** #
	id <- prepare_df_atc5(df_tot, atc_tokeep)

	# ******************************************************************** #
	# Prepare covariates to include in table one #
	# ******************************************************************** #

	#Socio-demographic variables and medical history
	explanatory1 = c("year_diag",
				"age","age_3cl",
				"depriv_index","depriv_index_quintile",
				"screening_before_diag",
				"gp_visit_4cl","gyneco_visit_4cl",
				"comedic_count_cl",
			    "charlson_indx_age")
	#Comorbid conditions
	df_comor <- id %>% select(starts_with("comor2"))
	explanatory2 =c(colnames(df_comor))
	#Other comedications
	df_comed <- id %>% select(matches("^[A-Z]{1}[0-9]{2}_bin$"),
							  matches(paste0("^",substr(atc_tokeep,1,3))))
	col_comed <- colnames(df_comed)
	explanatory3 =c(colnames(df_comed))
	#Dependent variable : medication intake
	dependent = "med_bin"
	df <- id %>% select(all_of(explanatory1),all_of(explanatory2),all_of(explanatory3),dependent)
	colnames_pretty <- data.frame(old_colnames = c(explanatory1,explanatory2,explanatory3)) %>%
		left_join(pretty_colnames)

	#Put label to columns
	for (col in setdiff(colnames(df),"med_bin")){
		print(col)
		if(!grepl("_bin|comor",col)){
			col_pretty <- as.character(colnames_pretty$pretty_colnames[colnames_pretty$old_colnames == col])
			col_pretty <- gsub("'"," ",col_pretty)
		}else if(grepl("comor2",col)){
			comor <- gsub("comor2_","",col)
			comor <- gsub("_"," ",comor)
			col_pretty <- comor
		}else{
			atc <- gsub("_bin","",col)
			col_pretty <- as.character(ref$PHA_ATC_LIB[ref$PHA_ATC_COD == atc])
			col_pretty <- paste0(col_pretty, " (",atc, ")")
		}
		print(col_pretty)
		text_to_parse <- paste0(
			"df <- df %>% mutate(",
			col,
			" = ff_label(",
			col,
			",'",
			col_pretty,
			"'))"
		)
		eval(parse(text = text_to_parse))
	}

	# ******************************************************************** #
	# Create table one (3 tables for explanatory 1, 2 and 3)
	# ******************************************************************** #
	t1 <- df %>%
		mutate(med_bin = ff_label(med_bin,name_atc_tokeep)) %>%
		summary_factorlist(dependent,explanatory1, p = TRUE, add_dependent_label = TRUE, dependent_label_prefix = "")
	t2 <- df %>%
		mutate(med_bin = ff_label(med_bin,name_atc_tokeep)) %>%
		summary_factorlist(dependent,explanatory2, p = TRUE, add_dependent_label = TRUE, dependent_label_prefix = "")
	t3 <- df %>%
		mutate(med_bin = ff_label(med_bin,name_atc_tokeep)) %>%
		summary_factorlist(dependent,explanatory3, p = TRUE, add_dependent_label = TRUE, dependent_label_prefix = "")

	# ******************************************************************** #
	# Save Table One
	# ******************************************************************** #

	#xlsx format
	t <- rbind(t1,t2,t3)
	write_xlsx(t, paste0(folder_txt,atc_tokeep,"_table1.xlsx"))

	#pdf format / one page per group of variables
	pdf(file = paste0(folder_txt,atc_tokeep,"_table1.pdf"), width=9,height=29, onefile = T)
		print(ggtexttable(t1,rows = NULL, theme = ttheme("light")) %>%
			tab_add_title(text = paste0("Table 1-a : Socio-demographic variables according to ",atc_tokeep, " consumption"), face = "bold"))
		print(ggtexttable(t2,rows = NULL, theme = ttheme("light")) %>%
			tab_add_title(text = paste0("Table 1-b : Comorbidites according to ",atc_tokeep, " consumption"), face = "bold"))
		print(ggtexttable(t3,rows = NULL, theme = ttheme("light")) %>%
			tab_add_title(text = paste0("Table 1-c : Other comedications according to ",atc_tokeep, " consumption"), face = "bold"))
	dev.off()

}

endsubmit;
run;

/* ******************************************************************** */
/* Step  3 : compute IPTW weights */
/* ******************************************************************** */

OPTIONS SET=R_HOME "C:\Program Files\R\R-3.6.3";
proc R;
export data = fresh.fresh_preprocessed_red R=df;  /* Database with all variables */
export data = fresh.top_atc5_snds R = list_comed;  /* List of medications to test */
submit;

#Libraries
library(tidyverse)
source("v10_util.R")
list_comed <- list_comed
df_tot <- df
#We will store all IPT weights in weight_atc5
weight_atc5 <- df %>% select(cletri)

#Loop on all medications
for (i in 1:nrow(list_comed)){

	#ATC code and label of the molecules
	atc_tokeep <- as.character(list_comed$PHA_ATC_C07[i]) <- atc
	name_atc_tokeep = list_comed$PHA_ATC_LIB[i] <- atc_name
	print(atc_tokeep)
	print(name_atc_tokeep)

	# ******************************************************************** #
	# Database preprocessing for the specific medication #
	# ******************************************************************** #
	id <- prepare_df_atc5(df_tot, atc_tokeep)

      # ******************************************************************** #
	  # Names of variables for comorbid conditions and other comedications
	  # The list of variables for other comedications may vary by comedication.
      # ******************************************************************** #

	#Comorbid conditions
	df_comor <- id %>% select(starts_with("comor2"))
	#Other comedications
	df_comed <- id %>% select(matches("^[A-Z]{1}[0-9]{2}_bin$"),
							  matches(paste0("^",substr(atc_tokeep,1,3))))
	col_comed <- colnames(df_comed)

    # ******************************************************************** #
	# Compute weights
    # ******************************************************************** #
	id <- compute_iptw_weights(id, df_comor, col_comed)

	# ******************************************************************** #
	# Store weights in database
    # ******************************************************************** #
	text_to_parse <- paste0(
		"weight_atc5$",
		atc_tokeep,
		"_iptw",
		"<- id$IPW_logistic_stab_trunc_right"
	)
	eval(parse(text = text_to_parse))
	text_to_parse <- paste0(
		"weight_atc5$",
		atc_tokeep,
		"_ps",
		"<- id$p_logistic_weights"
	)
	eval(parse(text = text_to_parse))

}

endsubmit;
import data = fresh.weight_atc5_work R = weight_atc5;
run;

/* We add weight atc2_work to weight_atc : enables to prevent erasing all database when code bugs*/
proc datasets library=fresh ; delete weight_atc5;
proc append base=fresh.weight_atc5 data=fresh.weight_atc5_work;

/* ******************************************************************** */
/* Step  4 : Check adjustment validity */
/* ******************************************************************** */

OPTIONS SET=R_HOME "C:\Program Files\R\R-3.6.3";
proc R;
export data = fresh.fresh_preprocessed_red R=df;  /* Database with all variables */
export data = fresh.top_atc5_snds R = list_comed;  /* List of medications to test */
export data = fresh.weight_atc5 R = weight; /* IPT weights */
submit;

#Libraries
library(tidyverse)
source("v10_util.R")

list_comed <- list_comed
df_tot <- df
df_tokeep_tot <- NULL

for (i in 1:nrow(list_comed)){

	#ATC code and label of the molecules
	atc_tokeep <- as.character(list_comed$PHA_ATC_C07[i]) <- atc
	name_atc_tokeep = list_comed$PHA_ATC_LIB[i] <- atc_name
	print(atc_tokeep)
	print(name_atc_tokeep)

	# ******************************************************************** #
	# Database preprocessing for the specific medication #
	# ******************************************************************** #
	id <- prepare_df_atc5(df_tot, atc_tokeep)
	#Add IPT weights to the database
	weight2 <- weight %>% select(cletri, p_logistic_weights = paste0(atc_tokeep,"_ps"),
								IPW_logistic_stab_trunc_right = paste0(atc_tokeep,"_iptw"))
	id <- id %>% left_join( weight2)

	# ************************************************************************* #
	# Compute Standardized Mean Differences for adjustment quality
	# ************************************************************************ #

	#Comorbid conditions variables
	df_comor <- id %>% select(starts_with("comor2"))
	#Other comedications variables
	df_comed <- id %>% select(matches("^[A-Z]{1}[0-9]{2}_bin$"),
							  matches(paste0("^",substr(atc_tokeep,1,3))))
	col_comed <- colnames(df_comed)
	#Compute and return SMD
	adj_quality = compute_adjustment_quality(id, df_comor, col_comed, atc_tokeep, atc_level="atc5")

	 # ******************************************************************** #
	 # Save adjustment quality summary #
	 # ******************************************************************** #
	tab_max = adj_quality$tab_max

	df_tokeep <- data.frame(
		atc_cod = atc_tokeep,
		atc_lib = name_atc_tokeep,
		n_atc = length(which(id$med_bin =="Yes")),
		worst_std_mean_diff = tab_max$Diff.Adj,
		var_worst_std_mean_diff = tab_max$Variable,
		adjustment_valid = ifelse(abs(as.numeric(tab_max$Diff.Adj)) > 0.1, "No","Yes")
		)
	df_tokeep_tot <- rbind(df_tokeep_tot,df_tokeep)

}

endsubmit;
import R= df_tokeep_tot data = fresh.atc5_adjustment_quality_work ;
run;

/* We add weight atc2_adjustment_quality_work to atc2_adjustment_quality : enables to prevent erasing all database when code bugs*/
proc datasets library=fresh ; delete atc5_adjustment_quality;
proc append base=fresh.atc5_adjustment_quality data=fresh.atc5_adjustment_quality_work;

/* ******************************************************************** */
/* Step  6 : Kaplan-Meier survival curves and ATE estimation
   for Overall Survival
   For whole population and by subgroup */
/* ******************************************************************** */

/* ******************************************************************** */
/* Step  6a : Draw survival curves */
/* ******************************************************************** */

OPTIONS SET=R_HOME "C:\Program Files\R\R-3.6.3";
proc R;
export data = fresh.fresh_preprocessed_red R=df;  /* Database with all variables */
export data = fresh.top_atc5_snds R = list_comed;  /* List of medications to test */
export data = fresh.weight_atc5 R = weight; /* IPT weights */
submit;

#Libraries
library(tidyverse)
source("v10_util.R")

list_comed <- list_comed
df_tot <- df
df_tokeep_tot <- NULL

#Loop on all medication
for (i in 1:nrow(list_comed)){

	#ATC code and label of the molecules
	atc_tokeep <- as.character(list_comed$atc_cod[i]) <- atc
	name_atc_tokeep = list_comed$atc_lib[i] <- atc_name
	print(atc_tokeep)
	print(name_atc_tokeep)

	# ******************************************************************** #
	# Database preprocessing for the specific medication #
	# ******************************************************************** #
	id <- prepare_df_atc5(df_tot, atc_tokeep)
	#Add IPT weights to the database
	weight2 <- weight %>% select(cletri, p_logistic_weights = paste0(atc_tokeep,"_ps"),
								IPW_logistic_stab_trunc_right = paste0(atc_tokeep,"_iptw"))
	id <- id %>% left_join( weight2)

	 # ******************************************************************** #
	 # Get the color for the medication (depends on ATC level 1 code)  #
	 # ******************************************************************** #
	color_palette <- c("#0DB14B","#CC004C","#CC004C","#F37021", "#FCB711","#0089D0","#6460AA","#006D77")
	names(color_palette) <- c("A","B","C","H","M","N","R","Others")
	atc_color <- as.character(color_palette[ifelse(substr(atc_tokeep,1,1) %in% names(color_palette), substr(atc_tokeep,1,1),"Others")])
	print(atc_color)

	 # ******************************************************************** #
	 # Create plot legend #
	 # ******************************************************************** #
	my_hist <- ggplot(data=id, aes(var,fill= var)) +
		geom_bar() +
		scale_fill_manual(values = c(atc_color,"#868686")) +
		labs(fill = name_atc_tokeep) +
		guides(fill = guide_legend(nrow = 1, title.position = "left"))
	legend <- cowplot::get_legend(my_hist)

	  # ******************************************************************** #
	  # Plots adjustment quality #
      # ******************************************************************** #

	#Density for weights
	plot_weights <- plot_ps(id,name_atc_tokeep,atc_color,atc_tokeep,atc_level=  "atc5")

	#Love plots for adjustment quality
	df_comor <- id %>% select(starts_with("comor2"))
	df_comed <- id %>% select(matches("^[A-Z]{1}[0-9]{2}_bin$"),
							  matches(paste0("^",substr(atc_tokeep,1,3))))
	col_comed <- colnames(df_comed)
	adj_quality = compute_adjustment_quality(id, df_comor, col_comed,atc_tokeep, atc_level = "atc5")

	#Adjustment quality check passed?
	tab_max = adj_quality$tab_max
	if(abs(as.numeric(tab_max$Diff.Adj))<= 0.1){
		decision_txt = "Adjustment quality check passed"
		decision_plot <- textGrob(decision_txt, gp = gpar(col = "darkgreen"))
	}else{
		decision_txt = "Adjustment quality check failed"
		decision_plot <- textGrob(decision_txt, gp = gpar(col = "red"))
	}

	#Merge plots for adjustment quality
	row6 <- plot_grid(plot_weights,adj_quality$love_plot2 + theme(legend.position = "none"),adj_quality$love_plot3,
			  ncol = 3 ,nrow = 1, labels = c('D','E','F'), rel_widths= c(3,3,3.5))
	row7 <- plot_grid(adj_quality$love_plot1 + theme(legend.position = "none"),
				adj_quality$love_plot4 + theme(legend.position = "none"),
				adj_quality$love_plot5,
			  ncol = 3 ,nrow = 1, labels = c('G','H','I'), rel_widths= c(3.5,2.5,3.5))

	  # ******************************************************************** #
	  # Kaplan-Meier survival curves
      # ******************************************************************** #

	#Whole population, chemotherapy yes, chemotherapy no
	data_ct <- id %>% filter(ct == "Yes")
	data_no_ct <- id %>% filter(ct == "No")
	row1 <- plot_survival_wp_ct_no_ct(id,atc,atc_name,atc_color,data_ct,data_no_ct,'J','K',atc_level = "atc5", estimate = "ATE")
	print(row1)

	#Endocrine therapy yes, endocrine therapy no
	data_ht <- id %>% filter(ht == "Yes")
	data_no_ht <- id %>% filter(ht == "No")
	plot_survival_endocrine_therapy(id,atc,atc_name,atc_color,data_ht,data_no_ht,atc_level ="atc5",estimate  ="ATE" )

	#Node positive, node negative
	data_nplus <- id %>% filter(pnuicc_2cl == "Node-positive")
	data_nmoins <- id %>% filter(pnuicc_2cl == "Node-negative")
	plot_survival_node_status(id,atc,atc_name,atc_color,data_nplus,
																data_nmoins,atc_level ="atc5",estimate  ="ATE")

	#By subtypes
	id_subtype <- id
	data_luminal <- id_subtype %>% filter(subtype == "luminal")
	data_her2 <- id_subtype %>% filter(subtype == "HER2+")
	data_tnbc <- id_subtype %>% filter(subtype == "TNBC")
	row234_col1 <- plot_survival_subtype(id,id_subtype,atc,atc_name,atc_color,
																data_luminal, data_her2,data_tnbc,atc_level = "atc5", estimate = "ATE")

     #By subtype in patients with chemotherapy
	data_ct_luminal <- id_subtype %>% filter(subtype == "luminal" & ct=="Yes")
	data_ct_her2 <- id_subtype %>% filter(subtype == "HER2+" & ct=="Yes")
	data_ct_tnbc <- id_subtype %>% filter(subtype == "TNBC" & ct=="Yes")
	row234_col2 <- plot_survival_subtype_ct(id,id_subtype,atc,atc_name,atc_color,
																data_ct_luminal, data_ct_her2,data_ct_tnbc,atc_level = "atc5", estimate = "ATE")

	#By subtype in patients without chemotherapy
	data_no_ct_luminal <- id_subtype %>% filter(subtype == "luminal" & ct=="No")
	data_no_ct_her2 <- id_subtype %>% filter(subtype == "HER2+" & ct=="No")
	data_no_ct_undefined <- id_subtype %>% filter(subtype == "undefined" & ct=="No")
	row234_col3 <- plot_survival_subtype_no_ct(id,id_subtype,atc,atc_name,atc_color,data_no_ct_luminal, data_no_ct_her2,data_no_ct_undefined,
																atc_level = "atc5",estimate = "ATE")

	 #By age classess
	 data_young <- id %>% filter(age_3cl == "<50")
	 data_middle <- id %>% filter(age_3cl == "50-75")
	 data_old <- id %>% filter(age_3cl == ">75")
	 row5 <- plot_survival_age(id,atc,atc_name,atc_color,data_young, data_middle,
															 data_old,atc_level = "atc5", estimate = "ATE")

	  # ******************************************************************** #
	  # Merge survival curves for whole population, chemotherapy, no chemotherapy
	  # subtype, subtype with chemotherapy, subtype without chemotherapy
	  # And save as pdf (in two pages)
      # ******************************************************************** #
		row234 <- plot_grid(row234_col1,row234_col2,row234_col3,
							ncol = 3 ,nrow = 1, labels = c('L','M',''))
		plot_total2 <- plot_grid(row6, row7, legend, as_ggplot(decision_plot), row1,
							ncol = 1 ,nrow = 5, rel_heights = c(1,1,0.1,0.1,1.3))
		plot_total3 <- plot_grid(row234, legend,
						ncol = 1 ,nrow = 2, rel_heights = c(3,0.2))

		pdf(file = paste0(folder_txt,atc,"_page23.pdf"), width=15.7,height=24, onefile = T)
			print(plot_total2)
			print(plot_total3)
		dev.off()

}

endsubmit;
run;

/* ******************************************************************** */
/* Step  6b : Compute Cox Hazard Ratio
# For whole population and by subgroup */
/* ********************************************************************* */

OPTIONS SET=R_HOME "C:\Program Files\R\R-3.6.3";
proc R;
export data = fresh.fresh_preprocessed_red R=df;  /* Database with all variables */
export data = fresh.top_atc5_snds R = list_comed;  /* List of medications to test */
export data = fresh.weight_atc5 R = weight; /* IPT weights */
submit;

#Libraries
library(tidyverse)
source("v10_util.R")

list_comed <- list_comed
df_tot <- df
df_tokeep_tot <- NULL

#Loop on all medications
for (i in 1:nrow(list_comed)){

	#ATC code and label of the molecules
	atc_tokeep <- as.character(list_comed$atc_cod[i]) <- atc
	name_atc_tokeep = list_comed$atc_lib[i] <- atc_name
	print(atc_tokeep)
	print(name_atc_tokeep)
	df_tokeep <- NULL

	# ******************************************************************** #
	# Database preprocessing for the specific medication #
	# ******************************************************************** #
	id <- prepare_df_atc5(df_tot, atc_tokeep)
	#Add IPT weights to the database
	weight2 <- weight %>% select(cletri, p_logistic_weights = paste0(atc_tokeep,"_ps"),
								IPW_logistic_stab_trunc_right = paste0(atc_tokeep,"_iptw"))
	id <- id %>% left_join( weight2)
	id$var <- factor(id$var,levels=c("No","Yes"))

	# ******************************************************************** #
	# Estimates for Whole Population (WP)  #
	# ******************************************************************** #

	info_tokeep_km <- data.frame("atc_cod" = atc_tokeep  ,
				"atc_lib"  = name_atc_tokeep)
	df_tokeep <- info_tokeep_km
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id,
					 weights = id$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_wp_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_wp_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_wp_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_wp_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_wp_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_wp_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	# ******************************************************************** #
	# Estimates by subgroup : endocrine therapy #
	# ******************************************************************** #

	#Compute interaction p-value
	mod1 <- coxph(Surv(delay_os, status_vital) ~  ht + var + ht*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	mod2 <- coxph(Surv(delay_os, status_vital) ~ ht + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	tmp_anova <- anova(mod2,mod1,test = "Chisq")
	info_tokeep_inter <- data.frame( "cox_p_interaction_med_ht" = tmp_anova[,"P(>|Chi|)"][2])
	df_tokeep <- cbind(df_tokeep,info_tokeep_inter)

	################
	# HT : yes
	################
	id_ht <- id %>% filter(ht == "Yes")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_ht,
					 weights = id_ht$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_ht_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_ht_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_ht_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_ht_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_ht_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_ht_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# HT : no
	################
	id_no_ht <- id %>% filter(ht == "No")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_no_ht,
					 weights = id_no_ht$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_no_ht_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_no_ht_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_no_ht_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_no_ht_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_no_ht_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_no_ht_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	# ******************************************************************** #
	# Estimates by subgroup : chemotherapy status #
	# ******************************************************************** #

	#Compute interaction p-value
	mod1 <- coxph(Surv(delay_os, status_vital) ~  ct + var + ct*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	mod2 <- coxph(Surv(delay_os, status_vital) ~ ct + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	tmp_anova <- anova(mod2,mod1,test = "Chisq")
	info_tokeep_inter <- data.frame( "cox_p_interaction_med_ct" = tmp_anova[,"P(>|Chi|)"][2])
	df_tokeep <- cbind(df_tokeep,info_tokeep_inter)

	################
	# CT : yes
	################
	id_ct <- id %>% filter(ct == "Yes")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_ct,
					 weights = id_ct$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_ct_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_ct_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_ct_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_ct_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_ct_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_ct_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# CT : no
	################
	id_no_ct <- id %>% filter(ct == "No")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_no_ct,
					 weights = id_no_ct$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_no_ct_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_no_ct_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_no_ct_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_no_ct_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_no_ct_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_no_ct_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	# ******************************************************************** #
	# Estimates by subgroup : nodal status #
	# ******************************************************************** #

	#Compute interaction p-value
	mod1 <- coxph(Surv(delay_os, status_vital) ~  pnuicc_2cl + var + pnuicc_2cl*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	mod2 <- coxph(Surv(delay_os, status_vital) ~ pnuicc_2cl + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	tmp_anova <- anova(mod2,mod1,test = "Chisq")
	info_tokeep_inter <- data.frame( "cox_p_interaction_med_node" = tmp_anova[,"P(>|Chi|)"][2])
	df_tokeep <- cbind(df_tokeep,info_tokeep_inter)

	################
	# Node-positive
	################
	id_nplus <- id %>% filter(pnuicc_2cl == "Node-positive")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_nplus,
					 weights = id_nplus$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_nplus_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_nplus_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_nplus_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_nplus_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_nplus_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_nplus_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# Node-negative
	################
	id_nmoins <- id %>% filter(pnuicc_2cl == "Node-negative")
	dim(id_nmoins)
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_nmoins,
					 weights = id_nmoins$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_nmoins_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_nmoins_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_nmoins_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_nmoins_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_nmoins_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_nmoins_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	# ******************************************************************** #
	# Estimates by subgroup : subtype #
	# ******************************************************************** #

	#Compute interaction p-value
	mod1 <- coxph(Surv(delay_os, status_vital) ~  subtype + var + subtype*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	mod2 <- coxph(Surv(delay_os, status_vital) ~ subtype + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	tmp_anova <- anova(mod2,mod1,test = "Chisq")
	info_tokeep_inter <- data.frame("cox_p_interaction_med_subtype" = tmp_anova[,"P(>|Chi|)"][2])
	df_tokeep <- cbind(df_tokeep,info_tokeep_inter)

	################
	# Luminal
	################
	id_luminal <- id %>% filter(subtype == "luminal")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_luminal,
					 weights = id_luminal$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_luminal_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_luminal_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_luminal_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_luminal_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_luminal_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_luminal_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# TNBC
	################
	id_tnbc <- id %>% filter(subtype == "TNBC")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_tnbc,
					 weights = id_tnbc$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph = tryCatch(expr = {cox.zph(fit = res.cox_weighted)$table[1,"p"]},
	error = function(e){print("error"); NA})
	if(is.na(sum$coefficients[1,"exp(coef)"])){
			info_tokeep_cox <- data.frame( "cox_coeff_tnbc_iptw" = NA,
				"cox_coeff_low_tnbc_iptw" = NA,
				"cox_coeff_high_tnbc_iptw" = NA,
				"cox_p_wald_tnbc_iptw" = NA ,
				"cox_p_robust_tnbc_iptw" =NA,
				"cox_zph_p_tnbc_iptw" = test_zph)
	}else{
			info_tokeep_cox <- data.frame( "cox_coeff_tnbc_iptw" =  sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_tnbc_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_tnbc_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_tnbc_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_tnbc_iptw" =sum$robscore[["pvalue"]],
				"cox_zph_p_tnbc_iptw" = test_zph)
	}
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# HER2+
	################
	id_her2 <- id %>% filter(subtype == "HER2+")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_her2,
					 weights = id_her2$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph = tryCatch(expr = {cox.zph(fit = res.cox_weighted)$table[1,"p"]},
	error = function(e){print("error"); NA})
	if(is.na(sum$coefficients[1,"exp(coef)"])){
			info_tokeep_cox <- data.frame( "cox_coeff_her2_iptw" = NA,
				"cox_coeff_low_her2_iptw" = NA,
				"cox_coeff_high_her2_iptw" = NA,
				"cox_p_wald_her2_iptw" = NA ,
				"cox_p_robust_her2_iptw" =NA,
				"cox_zph_p_her2_iptw" = test_zph)
	}else{
			info_tokeep_cox <- data.frame( "cox_coeff_her2_iptw" =  sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_her2_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_her2_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_her2_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_her2_iptw" =sum$robscore[["pvalue"]],
				"cox_zph_p_her2_iptw" = test_zph)
	}
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# Undefined
	################
	id_undefined <- id %>% filter(subtype == "undefined")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_undefined,
					 weights = id_undefined$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_undefined_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_undefined_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_undefined_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_undefined_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_undefined_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_undefined_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	# ******************************************************************** #
	# Estimates by subgroup : age #
	# ******************************************************************** #

	#Compute interaction p-value
	mod1 <- coxph(Surv(delay_os, status_vital) ~  age_3cl + var + age_3cl*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	mod2 <- coxph(Surv(delay_os, status_vital) ~ age_3cl + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	tmp_anova <- anova(mod2,mod1,test = "Chisq")
	info_tokeep_inter <- data.frame("cox_p_interaction_med_age" = tmp_anova[,"P(>|Chi|)"][2])
	df_tokeep <- cbind(df_tokeep,info_tokeep_inter)

	################
	# Young < 50
	################
	id_young <- id %>% filter(age_3cl == "<50")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_young,
					 weights = id_young$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_young_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_young_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_young_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_young_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_young_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_young_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# Intermediate : 50 - 75
	################
	id_middle <- id %>% filter(age_3cl == "50-75")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_middle,
					 weights = id_middle$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph = tryCatch(expr = {cox.zph(fit = res.cox_weighted)$table[1,"p"]},
	error = function(e){print("error"); NA})
	if(is.na(sum$coefficients[1,"exp(coef)"])){
			info_tokeep_cox <- data.frame( "cox_coeff_middle_iptw" = NA,
				"cox_coeff_low_middle_iptw" = NA,
				"cox_coeff_high_middle_iptw" = NA,
				"cox_p_wald_middle_iptw" = NA ,
				"cox_p_robust_middle_iptw" =NA,
				"cox_zph_p_middle_iptw" = test_zph)
	}else{
			info_tokeep_cox <- data.frame( "cox_coeff_middle_iptw" =  sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_middle_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_middle_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_middle_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_middle_iptw" =sum$robscore[["pvalue"]],
				"cox_zph_p_middle_iptw" = test_zph)
	}
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# Advanced : > 75
	################
	id_old <- id %>% filter(age_3cl == ">75")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_old,
					 weights = id_old$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph = tryCatch(expr = {cox.zph(fit = res.cox_weighted)$table[1,"p"]},
	error = function(e){print("error"); NA})
	if(is.na(sum$coefficients[1,"exp(coef)"])){
			info_tokeep_cox <- data.frame( "cox_coeff_old_iptw" = NA,
				"cox_coeff_low_old_iptw" = NA,
				"cox_coeff_high_old_iptw" = NA,
				"cox_p_wald_old_iptw" = NA ,
				"cox_p_robust_old_iptw" =NA,
				"cox_zph_p_old_iptw" = test_zph)
	}else{
			info_tokeep_cox <- data.frame( "cox_coeff_old_iptw" =  sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_old_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_old_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_old_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_old_iptw" =sum$robscore[["pvalue"]],
				"cox_zph_p_old_iptw" = test_zph)
	}
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# Store all results for this medication
	################
	df_tokeep_tot = rbind(df_tokeep_tot,df_tokeep)

}

endsubmit;
import data = fresh.atc5_ate_cox_work R = df_tokeep_tot;
run;

/* We add weight atc2cox_work to atc2cox : enables to prevent erasing all database when code bugs*/
proc datasets library=fresh ; delete atc5_ate_cox;
proc append base=fresh.atc5_ate_cox data=fresh.atc5_ate_cox_work;

/* Save ATE estimate  */
proc export data = fresh.atc5_ate_cox
	 outfile = "atc5_ate_cox.xlsx"
	dbms=xlsx
	replace;
run;

/* Save adjustment quality  */
proc export data = fresh.atc5_adjustment_quality
	 outfile = "atc5_adjusment_quality.xlsx"
	dbms=xlsx
	replace;
run;


/* ******************************************************************** */
/* Step  7 : Kaplan-Meier survival curves and ATE estimation
   for Disease free survival
   For whole population and by subgroup */
/* ******************************************************************** */

/* ******************************************************************** */
/* Step  7a : Draw survival curves */
/* ******************************************************************** */

OPTIONS SET=R_HOME "C:\Program Files\R\R-3.6.3";
proc R;
export data = fresh.dfs_preprocessed_red R=df;  /* Database with all variables */
export data = fresh.top_atc5_snds R = list_comed;  /* List of medications to test */
export data = fresh.weight_atc5 R = weight; /* IPT weights */
submit;

#Libraries
library(tidyverse)
source("v10_util.R")

list_comed <- list_comed
df_tot <- df
df_tokeep_tot <- NULL

#Loop on all medication
for (i in 1:nrow(list_comed)){

	#ATC code and label of the molecules
	atc_tokeep <- as.character(list_comed$atc_cod[i]) <- atc
	name_atc_tokeep = list_comed$atc_lib[i] <- atc_name
	print(atc_tokeep)
	print(name_atc_tokeep)

	# ******************************************************************** #
	# Database preprocessing for the specific medication #
	# ******************************************************************** #
	id <- prepare_df_atc5(df_tot, atc_tokeep)
	#Add IPT weights to the database
	weight2 <- weight %>% select(cletri, p_logistic_weights = paste0(atc_tokeep,"_ps"),
								IPW_logistic_stab_trunc_right = paste0(atc_tokeep,"_iptw"))
	id <- id %>% left_join( weight2)

	 # ******************************************************************** #
	 # Get the color for the medication (depends on ATC level 1 code)  #
	 # ******************************************************************** #
	color_palette <- c("#0DB14B","#CC004C","#CC004C","#F37021", "#FCB711","#0089D0","#6460AA","#006D77")
	names(color_palette) <- c("A","B","C","H","M","N","R","Others")
	atc_color <- as.character(color_palette[ifelse(substr(atc_tokeep,1,1) %in% names(color_palette), substr(atc_tokeep,1,1),"Others")])
	print(atc_color)

	 # ******************************************************************** #
	 # Create plot legend #
	 # ******************************************************************** #
	my_hist <- ggplot(data=id, aes(var,fill= var)) +
		geom_bar() +
		scale_fill_manual(values = c(atc_color,"#868686")) +
		labs(fill = name_atc_tokeep) +
		guides(fill = guide_legend(nrow = 1, title.position = "left"))
	legend <- cowplot::get_legend(my_hist)

	  # ******************************************************************** #
	  # Plots adjustment quality #
      # ******************************************************************** #

	#Density for weights
	plot_weights <- plot_ps(id,name_atc_tokeep,atc_color,atc_tokeep,atc_level=  "atc5",estimate = "ATE_rfs")

	#Love plots for adjustment quality
	df_comor <- id %>% select(starts_with("comor2"))
	df_comed <- id %>% select(matches("^[A-Z]{1}[0-9]{2}_bin$"),
							  matches(paste0("^",substr(atc_tokeep,1,3))))
	col_comed <- colnames(df_comed)
	adj_quality = compute_adjustment_quality(id, df_comor, col_comed,atc_tokeep, atc_level = "atc5",estimate = "ATE_rfs")

	#Adjustment quality check passed?
	tab_max = adj_quality$tab_max
	if(abs(as.numeric(tab_max$Diff.Adj))<= 0.1){
		decision_txt = "Adjustment quality check passed"
		decision_plot <- textGrob(decision_txt, gp = gpar(col = "darkgreen"))
	}else{
		decision_txt = "Adjustment quality check failed"
		decision_plot <- textGrob(decision_txt, gp = gpar(col = "red"))
	}

	#Merge plots for adjustment quality
	row6 <- plot_grid(plot_weights,adj_quality$love_plot2 + theme(legend.position = "none"),adj_quality$love_plot3,
			  ncol = 3 ,nrow = 1, labels = c('D','E','F'), rel_widths= c(3,3,3.5))
	row7 <- plot_grid(adj_quality$love_plot1 + theme(legend.position = "none"),
				adj_quality$love_plot4 + theme(legend.position = "none"),
				adj_quality$love_plot5,
			  ncol = 3 ,nrow = 1, labels = c('G','H','I'), rel_widths= c(3.5,2.5,3.5))

	  # ******************************************************************** #
	  # Kaplan-Meier survival curves
      # ******************************************************************** #

	#Whole population, chemotherapy yes, chemotherapy no
	data_ct <- id %>% filter(ct == "Yes")
	data_no_ct <- id %>% filter(ct == "No")
	row1 <- plot_survival_wp_ct_no_ct(id,atc,atc_name,atc_color,data_ct,data_no_ct,'J','K',atc_level = "atc5",
				estimate = "ATE_rfs", ylab_txt = "Disease-free survival probability")
	print(row1)

	#Endocrine therapy yes, endocrine therapy no
	data_ht <- id %>% filter(ht == "Yes")
	data_no_ht <- id %>% filter(ht == "No")
	plot_survival_endocrine_therapy(id,atc,atc_name,atc_color,data_ht,data_no_ht,atc_level ="atc5",estimate  ="ATE_rfs",
												ylab_txt = "Disease-free survival probability")

	#Node positive, node negative
	data_nplus <- id %>% filter(pnuicc_2cl == "Node-positive")
	data_nmoins <- id %>% filter(pnuicc_2cl == "Node-negative")
	plot_survival_node_status(id,atc,atc_name,atc_color,data_nplus,
																data_nmoins,atc_level ="atc5",estimate  ="ATE_rfs",
																ylab_txt = "Disease-free survival probability")

	#By subtypes
	id_subtype <- id
	data_luminal <- id_subtype %>% filter(subtype == "luminal")
	data_her2 <- id_subtype %>% filter(subtype == "HER2+")
	data_tnbc <- id_subtype %>% filter(subtype == "TNBC")
	row234_col1 <- plot_survival_subtype(id,id_subtype,atc,atc_name,atc_color,
																data_luminal, data_her2,data_tnbc,atc_level = "atc5",
																estimate = "ATE_rfs", ylab_txt = "Disease-free survival probability")

     #By subtype in patients with chemotherapy
	data_ct_luminal <- id_subtype %>% filter(subtype == "luminal" & ct=="Yes")
	data_ct_her2 <- id_subtype %>% filter(subtype == "HER2+" & ct=="Yes")
	data_ct_tnbc <- id_subtype %>% filter(subtype == "TNBC" & ct=="Yes")
	row234_col2 <- plot_survival_subtype_ct(id,id_subtype,atc,atc_name,atc_color,
																data_ct_luminal, data_ct_her2,data_ct_tnbc,atc_level = "atc5", estimate = "ATE_rfs", ylab_txt = "Disease-free survival probability")

	#By subtype in patients without chemotherapy
	data_no_ct_luminal <- id_subtype %>% filter(subtype == "luminal" & ct=="No")
	data_no_ct_her2 <- id_subtype %>% filter(subtype == "HER2+" & ct=="No")
	data_no_ct_undefined <- id_subtype %>% filter(subtype == "undefined" & ct=="No")
	row234_col3 <- plot_survival_subtype_no_ct(id,id_subtype,atc,atc_name,atc_color,data_no_ct_luminal, data_no_ct_her2,data_no_ct_undefined,
																atc_level = "atc5",estimate = "ATE_rfs", ylab_txt = "Disease-free survival probability")

	 #By age classess
	 data_young <- id %>% filter(age_3cl == "<50")
	 data_middle <- id %>% filter(age_3cl == "50-75")
	 data_old <- id %>% filter(age_3cl == ">75")
	 row5 <- plot_survival_age(id,atc,atc_name,atc_color,data_young, data_middle,
															 data_old,atc_level = "atc5", estimate = "ATE_rfs", ylab_txt = "Disease-free survival probability")

	  # ******************************************************************** #
	  # Merge survival curves for whole population, chemotherapy, no chemotherapy
	  # subtype, subtype with chemotherapy, subtype without chemotherapy
	  # And save as pdf (in two pages)
      # ******************************************************************** #
		row234 <- plot_grid(row234_col1,row234_col2,row234_col3,
							ncol = 3 ,nrow = 1, labels = c('L','M',''))
		plot_total2 <- plot_grid(row6, row7, legend, as_ggplot(decision_plot), row1,
							ncol = 1 ,nrow = 5, rel_heights = c(1,1,0.1,0.1,1.3))
		plot_total3 <- plot_grid(row234, legend,
						ncol = 1 ,nrow = 2, rel_heights = c(3,0.2))

		pdf(file = paste0(folder_txt,atc,"_page23.pdf"), width=15.7,height=24, onefile = T)
			print(plot_total2)
			print(plot_total3)
		dev.off()

}

endsubmit;
run;

/* ******************************************************************** */
/* Step  7b : Compute Cox Hazard Ratio
# For whole population and by subgroup */
/* ********************************************************************* */

OPTIONS SET=R_HOME "C:\Program Files\R\R-3.6.3";
proc R;
export data = fresh.dfs_preprocessed_red R=df;  /* Database with all variables */
export data = fresh.top_atc5_snds R = list_comed;  /* List of medications to test */
export data = fresh.weight_atc5 R = weight; /* IPT weights */
submit;

#Libraries
library(tidyverse)
source("v10_util.R")

list_comed <- list_comed
df_tot <- df
df_tokeep_tot <- NULL

#Loop on all medications
for (i in 1:nrow(list_comed)){

	#ATC code and label of the molecules
	atc_tokeep <- as.character(list_comed$atc_cod[i]) <- atc
	name_atc_tokeep = list_comed$atc_lib[i] <- atc_name
	print(atc_tokeep)
	print(name_atc_tokeep)
	df_tokeep <- NULL

	# ******************************************************************** #
	# Database preprocessing for the specific medication #
	# ******************************************************************** #
	id <- prepare_df_atc5(df_tot, atc_tokeep)
	#Add IPT weights to the database
	weight2 <- weight %>% select(cletri, p_logistic_weights = paste0(atc_tokeep,"_ps"),
								IPW_logistic_stab_trunc_right = paste0(atc_tokeep,"_iptw"))
	id <- id %>% left_join( weight2)
	id$var <- factor(id$var,levels=c("No","Yes"))

	# ******************************************************************** #
	# Estimates for Whole Population (WP)  #
	# ******************************************************************** #

	info_tokeep_km <- data.frame("atc_cod" = atc_tokeep  ,
				"atc_lib"  = name_atc_tokeep)
	df_tokeep <- info_tokeep_km
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id,
					 weights = id$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_wp_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_wp_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_wp_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_wp_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_wp_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_wp_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	# ******************************************************************** #
	# Estimates by subgroup : endocrine therapy #
	# ******************************************************************** #

	#Compute interaction p-value
	mod1 <- coxph(Surv(delay_os, status_vital) ~  ht + var + ht*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	mod2 <- coxph(Surv(delay_os, status_vital) ~ ht + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	tmp_anova <- anova(mod2,mod1,test = "Chisq")
	info_tokeep_inter <- data.frame( "cox_p_interaction_med_ht" = tmp_anova[,"P(>|Chi|)"][2])
	df_tokeep <- cbind(df_tokeep,info_tokeep_inter)

	################
	# HT : yes
	################

	id_ht <- id %>% filter(ht == "Yes")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_ht,
					 weights = id_ht$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_ht_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_ht_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_ht_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_ht_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_ht_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_ht_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# HT : no
	################
	id_no_ht <- id %>% filter(ht == "No")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_no_ht,
					 weights = id_no_ht$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_no_ht_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_no_ht_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_no_ht_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_no_ht_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_no_ht_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_no_ht_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	# ******************************************************************** #
	# Estimates by subgroup : chemotherapy status #
	# ******************************************************************** #

	#Compute interaction p-value
	mod1 <- coxph(Surv(delay_os, status_vital) ~  ct + var + ct*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	mod2 <- coxph(Surv(delay_os, status_vital) ~ ct + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	tmp_anova <- anova(mod2,mod1,test = "Chisq")
	info_tokeep_inter <- data.frame( "cox_p_interaction_med_ct" = tmp_anova[,"P(>|Chi|)"][2])
	df_tokeep <- cbind(df_tokeep,info_tokeep_inter)

	################
	# CT : yes
	################
	id_ct <- id %>% filter(ct == "Yes")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_ct,
					 weights = id_ct$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_ct_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_ct_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_ct_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_ct_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_ct_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_ct_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# CT : no
	################
	id_no_ct <- id %>% filter(ct == "No")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_no_ct,
					 weights = id_no_ct$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_no_ct_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_no_ct_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_no_ct_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_no_ct_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_no_ct_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_no_ct_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	# ******************************************************************** #
	# Estimates by subgroup : nodal status #
	# ******************************************************************** #

	#Compute interaction p-value
	mod1 <- coxph(Surv(delay_os, status_vital) ~  pnuicc_2cl + var + pnuicc_2cl*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	mod2 <- coxph(Surv(delay_os, status_vital) ~ pnuicc_2cl + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	tmp_anova <- anova(mod2,mod1,test = "Chisq")
	info_tokeep_inter <- data.frame( "cox_p_interaction_med_node" = tmp_anova[,"P(>|Chi|)"][2])
	df_tokeep <- cbind(df_tokeep,info_tokeep_inter)

	################
	# Node-positive
	################
	id_nplus <- id %>% filter(pnuicc_2cl == "Node-positive")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_nplus,
					 weights = id_nplus$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_nplus_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_nplus_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_nplus_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_nplus_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_nplus_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_nplus_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# Node-negative
	################
	id_nmoins <- id %>% filter(pnuicc_2cl == "Node-negative")
	dim(id_nmoins)
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_nmoins,
					 weights = id_nmoins$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_nmoins_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_nmoins_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_nmoins_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_nmoins_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_nmoins_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_nmoins_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	# ******************************************************************** #
	# Estimates by subgroup : subtype #
	# ******************************************************************** #

	#Compute interaction p-value
	mod1 <- coxph(Surv(delay_os, status_vital) ~  subtype + var + subtype*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	mod2 <- coxph(Surv(delay_os, status_vital) ~ subtype + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	tmp_anova <- anova(mod2,mod1,test = "Chisq")
	info_tokeep_inter <- data.frame("cox_p_interaction_med_subtype" = tmp_anova[,"P(>|Chi|)"][2])
	df_tokeep <- cbind(df_tokeep,info_tokeep_inter)

	################
	# Luminal
	################
	id_luminal <- id %>% filter(subtype == "luminal")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_luminal,
					 weights = id_luminal$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_luminal_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_luminal_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_luminal_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_luminal_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_luminal_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_luminal_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# TNBC
	################
	id_tnbc <- id %>% filter(subtype == "TNBC")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_tnbc,
					 weights = id_tnbc$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph = tryCatch(expr = {cox.zph(fit = res.cox_weighted)$table[1,"p"]},
	error = function(e){print("error"); NA})
	if(is.na(sum$coefficients[1,"exp(coef)"])){
			info_tokeep_cox <- data.frame( "cox_coeff_tnbc_iptw" = NA,
				"cox_coeff_low_tnbc_iptw" = NA,
				"cox_coeff_high_tnbc_iptw" = NA,
				"cox_p_wald_tnbc_iptw" = NA ,
				"cox_p_robust_tnbc_iptw" =NA,
				"cox_zph_p_tnbc_iptw" = test_zph)
	}else{
			info_tokeep_cox <- data.frame( "cox_coeff_tnbc_iptw" =  sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_tnbc_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_tnbc_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_tnbc_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_tnbc_iptw" =sum$robscore[["pvalue"]],
				"cox_zph_p_tnbc_iptw" = test_zph)
	}
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# HER2+
	################
	id_her2 <- id %>% filter(subtype == "HER2+")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_her2,
					 weights = id_her2$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph = tryCatch(expr = {cox.zph(fit = res.cox_weighted)$table[1,"p"]},
	error = function(e){print("error"); NA})
	if(is.na(sum$coefficients[1,"exp(coef)"])){
			info_tokeep_cox <- data.frame( "cox_coeff_her2_iptw" = NA,
				"cox_coeff_low_her2_iptw" = NA,
				"cox_coeff_high_her2_iptw" = NA,
				"cox_p_wald_her2_iptw" = NA ,
				"cox_p_robust_her2_iptw" =NA,
				"cox_zph_p_her2_iptw" = test_zph)
	}else{
			info_tokeep_cox <- data.frame( "cox_coeff_her2_iptw" =  sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_her2_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_her2_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_her2_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_her2_iptw" =sum$robscore[["pvalue"]],
				"cox_zph_p_her2_iptw" = test_zph)
	}
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# Undefined
	################
	id_undefined <- id %>% filter(subtype == "undefined")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_undefined,
					 weights = id_undefined$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_undefined_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_undefined_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_undefined_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_undefined_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_undefined_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_undefined_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	# ******************************************************************** #
	# Estimates by subgroup : age #
	# ******************************************************************** #

	#Compute interaction p-value
	mod1 <- coxph(Surv(delay_os, status_vital) ~  age_3cl + var + age_3cl*var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	mod2 <- coxph(Surv(delay_os, status_vital) ~ age_3cl + var, data = id, weights = id$IPW_logistic_stab_trunc_right, robust = TRUE)
	tmp_anova <- anova(mod2,mod1,test = "Chisq")
	info_tokeep_inter <- data.frame("cox_p_interaction_med_age" = tmp_anova[,"P(>|Chi|)"][2])
	df_tokeep <- cbind(df_tokeep,info_tokeep_inter)

	################
	# Young < 50
	################
	id_young <- id %>% filter(age_3cl == "<50")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_young,
					 weights = id_young$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph  = cox.zph(fit = res.cox_weighted)
	info_tokeep_cox <- data.frame( "cox_coeff_young_iptw" = sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_young_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_young_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_young_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_young_iptw" = sum$robscore[["pvalue"]],
				"cox_zph_p_young_iptw" = test_zph$table[1,"p"])
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# Intermediate : 50 - 75
	################
	id_middle <- id %>% filter(age_3cl == "50-75")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_middle,
					 weights = id_middle$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph = tryCatch(expr = {cox.zph(fit = res.cox_weighted)$table[1,"p"]},
	error = function(e){print("error"); NA})
	if(is.na(sum$coefficients[1,"exp(coef)"])){
			info_tokeep_cox <- data.frame( "cox_coeff_middle_iptw" = NA,
				"cox_coeff_low_middle_iptw" = NA,
				"cox_coeff_high_middle_iptw" = NA,
				"cox_p_wald_middle_iptw" = NA ,
				"cox_p_robust_middle_iptw" =NA,
				"cox_zph_p_middle_iptw" = test_zph)
	}else{
			info_tokeep_cox <- data.frame( "cox_coeff_middle_iptw" =  sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_middle_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_middle_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_middle_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_middle_iptw" =sum$robscore[["pvalue"]],
				"cox_zph_p_middle_iptw" = test_zph)
	}
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# Advanced : > 75
	################
	id_old <- id %>% filter(age_3cl == ">75")
	#Cox Hazard Ratio
	res.cox_weighted = coxph(Surv(delay_os,status_vital)  ~ var, data = id_old,
					 weights = id_old$IPW_logistic_stab_trunc_right,
					 robust = TRUE)
	sum <- summary(res.cox_weighted)
	test_zph = tryCatch(expr = {cox.zph(fit = res.cox_weighted)$table[1,"p"]},
	error = function(e){print("error"); NA})
	if(is.na(sum$coefficients[1,"exp(coef)"])){
			info_tokeep_cox <- data.frame( "cox_coeff_old_iptw" = NA,
				"cox_coeff_low_old_iptw" = NA,
				"cox_coeff_high_old_iptw" = NA,
				"cox_p_wald_old_iptw" = NA ,
				"cox_p_robust_old_iptw" =NA,
				"cox_zph_p_old_iptw" = test_zph)
	}else{
			info_tokeep_cox <- data.frame( "cox_coeff_old_iptw" =  sum$coefficients[1,"exp(coef)"],
				"cox_coeff_low_old_iptw" = sum$conf.int[1,"lower .95"],
				"cox_coeff_high_old_iptw" = sum$conf.int[1,"upper .95"],
				"cox_p_wald_old_iptw" = sum$waldtest[["pvalue"]] ,
				"cox_p_robust_old_iptw" =sum$robscore[["pvalue"]],
				"cox_zph_p_old_iptw" = test_zph)
	}
	df_tokeep <- cbind(df_tokeep,info_tokeep_cox)

	################
	# Store all results for this medication
	################
	df_tokeep_tot = rbind(df_tokeep_tot,df_tokeep)

}

endsubmit;
import data = fresh.atc5_aterfs_work R = df_tokeep_tot;
run;

/* We add weight atc5cox_work to atc5cox : enables to prevent erasing all database when code bugs*/
proc datasets library=fresh ; delete atc5_aterfs;
proc append base=fresh.atc5_aterfs data=fresh.atc5_aterfs_work;

/* Save ATE estimate */
proc export data = fresh.atc5_aterfs
	 outfile = "atc5_ate_cox.xlsx"
	dbms=xlsx
	replace;
run;
