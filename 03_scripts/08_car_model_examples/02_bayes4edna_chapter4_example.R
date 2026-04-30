## - Bayesian models for eDNA (2023 samples)


## created 5 July 2024
## by Molly M Kressler

########
## Load data  
########
## cleaned in pipeline documents

pacman::p_load(sf,tidyverse,dplyr,ggplot2, patchwork, cowplot,lubridate,flextable, rnaturalearth, lme4, modelsummary, readr, readxl, brms, bayesplot, tidybayes, modelr, ggdist, nimble, MCMCvis)

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/')

## qPCR results, all species, not standards, test assay positive nor negative test controls 
fieldsamples <- read.csv('EDNA/data_edna/qPCRresults/processedQPCRresults_cornwallspecies_june2024.csv')%>%
	filter(NTC.Amplified == 0)%>%
	dplyr::select(-NTC.Amplified)%>%
		as_tibble() # st_write converted the TRUE/FALSE to 0s and 1s. 0 = TRUE and 1 = FALSE. This removes any plate where the NTC amplified before the LOQ. 

########
## Data exploration 
########

	hist<- ggplot(data = fieldsamples, aes(x =log(copies.sampavg+1), fill = Target.Name))+
		geom_histogram(position = 'dodge', bins = 5)+
	  	scale_fill_manual(values = c('#477939', '#799ecb', '#85cb7c', '#003a78'), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus'))+
		theme_bw()+
		labs(y = NULL, x = 'DNA yield (log10)')

	ggsave(hist, file = 'EDNA/data_edna/brms/figures/hist_bySpecies.png', device = 'png', units = 'in', height = 6, width = 6, dpi = 800)


########
## BRMs models - normal distribution with log10 transformation 
########

	## prepare data 
	data <- fieldsamples %>% 
		mutate(taxa = as.character(case_when(Target.Name == 'Engraulis' ~ 'Fish', Target.Name == 'Scomber' ~ 'Fish',Target.Name == 'Prionace' ~ 'Shark', Target.Name == 'Alopias' ~ 'Shark')))%>%
		mutate(log10.copies = log(copies.sampavg+0.1),
			round.copies = as.numeric(round(copies.sampavg,0)),
			re.level = as.numeric(as.numeric(factor(eventID))),
			methodID = as.numeric(as.numeric(factor(method.type))),
			taxaID = as.numeric(as.numeric(factor(taxa))))%>%
		dplyr::select('eventID', 're.level', 'Sample.Name', 'Target.Name','taxa', 'method.type', 'methodID', 'taxaID','copies.techrepavg','log10.copies', 'round.copies')
	data # waterbottles = 2, metaprobes = 1,fish = 1, shark = 2

	######
	### Fit Models
	######		

	## fit1 = gaussian with log10 data
	fit1 <- brm(log10.copies ~ method.type * taxa + (1|eventID), 
			data = data, 
			family = gaussian(), 
			file = 'EDNA/data_edna/brms/model1_gaussian_log10data.RDS') # auto-saves the model object as an RDS file

	fit1 <- readRDS('EDNA/data_edna/brms/model1_gaussian_log10data.RDS')

	## fit2 = poisson - hates this. wont sample. 
	fit2 <- brm(round.copies ~ method.type * taxa + (1|eventID), 
			data = data, 
			family = poisson(link = 'log'), 
			file = 'EDNA/data_edna/brms/model2_Poisson_notransformation.RDS') # auto-saves the model object as an RDS file


########
## Nimble  models - poissons 
########
	# nimble model for poisson fit2 
		modelcode_fit2 <- nimbleCode({

			# priors - taxa #
			for(i in 1:K){
				m[i] ~ dnorm(0,.01) 
				t[i] ~ dnorm(0,.01) 
			}
			# priors for random intercept, informative with Half-Cauchy priors for sigma
			for(i in 1:E){
				e[i] ~ dnorm(0, tau.re)
			}
				num ~ dnorm(0, 0.0016)
				denom ~ dnorm(0,1)
				sigma.re <- abs(num/denom)
				tau.re <- 1/(sigma.re * sigma.re)
			
			# likelihood #
			for(i in 1:N){
				Y[i] ~ dpois(mu[i])
				log(mu[i]) <- (m[methodID[i]]*methodID[i]) * (t[taxaID[i]]*taxaID[i]) + e[re[i]]
				#log(mu[i]) <- m[methodID[i]] * t[taxaID[i]] + e[re[i]]

				# Discrepancy measures
				#Ynew[i] ~ dpois(mu[i]) #new data 
				#ExpY[i] <- mu[i]
				#VarY[i] <- mu[i]
				#PRes[i] <- (Y[i] - ExpY[i])/sqrt(VarY[i])
				#PResNew[i] <- (Ynew[i] - ExpY[i])/sqrt(VarY[i])
				#D[i] <- pow(PRes[i], 2)
				#Dnew[i] <- pow(PResNew[i], 2)
			}
		})
		
		### no random effect
		
		modelcode_fit2_noRE <- nimbleCode({
			# priors - taxa #
			for(i in 1:K){
				m[i] ~ dnorm(0, 0.01) 
				t[i] ~ dnorm(0, 0.01) 
			}
						
			# likelihood #
			for(i in 1:N){
				Y[i] ~ dpois(mu[i])
				log(mu[i]) <- (m[methodID[i]]*methodID[i]) * (t[taxaID[i]]*taxaID[i]) 
				#log(mu[i]) <- m[methodID[i]] * t[taxaID[i]] 

			}
		})


		### bits and bobs for compilation
		fit2.constants <- list(
				K = max(data$methodID),
				E = length(unique(data$re.level)), 
				N = nrow(data),
				methodID = data$methodID, 
				taxaID = data$taxaID,
				re = data$re.level)

		fit2.data <- list(Y = round(data$copies.techrepavg,0))

		fit2.init <- list(
			e = rnorm(15,0,1),
			m = rnorm(2, 0, 1),
			t = rnorm(2, 0, 1),
			tau.re = rgamma(1,1,1),
			num = rnorm(1, 0, 0.0016),
			denom = rnorm(1, 0, 1)
			)

		### compile
	    model_fit2<-nimbleModel(code=modelcode_fit2, name="model_fit2",data=fit2.data, constants = fit2.constants,inits=fit2.init) #define the model

	    	model_fit2$calculate() # if = NA, indicates missing or invalid initial values, and you have to fix the model until it is numeric.
	    	
	    	Cmf2<-compileNimble(model_fit2) 
	    	conf2 <- configureMCMC(model_fit2, monitors = c('tau.re', 'm', 't', 'e'), onlySLICE = FALSE)
		    MCMC_model_fit2 <- buildMCMC(conf2, na.rm = TRUE)
		    ccMCMC2 <- compileNimble(MCMC_model_fit2, project = model_fit2)
		    samples2 <- runMCMC(ccMCMC2, niter = 65000, nburnin = 2500, nchains = 3, samplesAsCodaMCMC = TRUE)

	    model_fit2_noRE<-nimbleModel(code=modelcode_fit2_noRE, name="model_fit2",data=fit2.data, constants = fit2.constants,inits=fit2.init) #define the model
	    	Cmf2_noRE<-compileNimble(model_fit2_noRE) 
	    	conf2_noRE <- configureMCMC(model_fit2_noRE, monitors = c( 'm', 't'), onlySLICE = FALSE)
		    MCMC_model_fit2_noRE <- buildMCMC(conf2_noRE, na.rm = TRUE)
		    ccMCMC2_noRE <- compileNimble(MCMC_model_fit2_noRE, project = model_fit2_noRE)
		    samples2_noRE <- runMCMC(ccMCMC2_noRE, niter = 85000, nburnin = 8500, nchains = 3, samplesAsCodaMCMC = TRUE)
	    
	    MCMCsummary(samples2)
	    MCMCsummary(samples2_noRE)

	######
	### Exploration
	######	
	get_variables(fit1) # raw list of variable names. b_Intercept is the global mean, and the r_eventID[level, Intercept] are offsets of that mean for each condition. 

	## flextable of summary
	fit1table <- modelsummary(fit1, fmt=3,estimate='estimate', statistic='conf.int',conf_level=0.95,output='flextable')%>%
		theme_zebra()%>%
		set_header_labels('Model 1' = 'Method & Target\n\ Taxa GLMM')%>%
	    align(align = 'center', part = 'all')%>%
	    font(fontname = 'Arial', part = 'all')%>%
	    fontsize(size = 10, part = 'all')%>%
	    autofit()

	    ## fit 2
	    	samples2
		    summary_fit2 <- MCMCsummary(samples2) %>%
		      tibble::rownames_to_column()%>%
		      rename_with(str_to_title)%>%
		      #rename(Parameter = Rowname)%>%
		      #rename('% of posterior with \n\ same sign as estimate' = 'P>0', Estimate = 'Mean','lower'='5%',upper='95%')%>%
		      #mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
		      #dplyr::select(-lower,-upper,-Sd)%>% 
		      flextable()%>%
		      theme_zebra()%>%
		      set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
		      align(align = 'center', part = 'all')%>%
		      font(fontname = 'Arial', part = 'all')%>%
		      color(color='black',part='all')%>%
		      fontsize(size = 10, part = 'all')%>%
		      autofit()
		    summary_fit2


	# Gather variable indices into a separate column (one row per draw, one column per parameter)
	## fit 1
		draws <- fit1 %>%
			gather_draws(`b_Intercept`, b_method.typewaterbottle, b_taxaShark,`b_method.typewaterbottle:taxaShark`, r_eventID[condition, ])%>%
			median_qi(.width = c(.95, .5)) # calculates the median and the 95% quantile interval for the parameter

		caterpillars <- draws %>%
			dplyr::filter(.variable != 'r_eventID') %>%
			ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper))+
			geom_pointinterval()+      
			geom_vline(xintercept=0,linetype=3)+
			theme_bw()
		caterpillars

		ggsave(caterpillars, file = 'EDNA/data_edna/brms/figures/caterpillars_model1_gaussian_log10data.png', device = 'png', unit = 'in', height = 4, dpi = 850)	

	## fit 2
		draws <- model_fit2 %>%
			gather_draws(`b_Intercept`, b_method.typewaterbottle, b_taxaShark,`b_method.typewaterbottle:taxaShark`, r_eventID[condition, ])%>%
			median_qi(.width = c(.95, .5)) # calculates the median and the 95% quantile interval for the parameter

		caterpillars <- draws %>%
			dplyr::filter(.variable != 'r_eventID') %>%
			ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper))+
			geom_pointinterval()+      
			geom_vline(xintercept=0,linetype=3)+
			theme_bw()
		caterpillars

		ggsave(caterpillars, file = 'EDNA/data_edna/brms/figures/caterpillars_model1_gaussian_log10data.png', device = 'png', unit = 'in', height = 4, dpi = 850)	



	# Add draws from posterior fit and (then) residuals, using epred_draws

	postfit <- data %>%
		dplyr::select('eventID', 'Sample.Name', 'taxa', 'method.type', 'log10.copies')%>%
		data_grid(taxa, method.type, eventID)%>%
		add_epred_draws(fit1)%>%
		ggplot(aes(x = .epred, y = interaction(taxa,method.type), fill = method.type))+
		stat_halfeye(alpha = 0.75)+
		scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
		geom_vline(xintercept=0,linetype=3)+
		labs(y = 'Interaction of Taxa & Method', x = 'Conditional Means')+
		guides(fill = 'none')+
		theme_bw()
	ggsave(postfit, file = 'EDNA/data_edna/brms/figures/postdraws_caterpillars_fixedeffects_model1_gaussian_log10data.png', device = 'png', unit = 'in', height = 5, dpi = 850)	


	## Diagnostics 

		## Overdispersion 
		n <- nrow(data)
			## fit1 
			res_fit1 <- resid(fit1)
			overdis_fit1 = sum(res_fit1^2)/(n-2)
			overdis_fit1 # 18.915, overdispersed (Gaussian on log10 copies)

		## Plots of residuals

			## fit 1
			res_fit1 <- as_tibble(resid(fit1)) # E1 in the book, page 58
			F1 <- as_tibble(predict(fit1))%>%rename(fitted = Estimate, fit.Q2.5 = Q2.5, fitQ97.5 = Q97.5, fitEst.Error = Est.Error)
			diag_fit1 <- bind_cols(F1, res_fit1)

			resVfit_fit1 <- ggplot(data = diag_fit1,aes(x= fitted, y = Estimate))+
				geom_point()+ 
				labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
				theme_bw()
			res_fit1_hist <- ggplot(data = res_fit1, aes(x = Estimate))+ 
				geom_histogram(binwidth = nrow(res_fit1)/1000,fill = 'black')+ 
				theme_bw()+
				labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
			res_fit1_qq <- ggplot(data=F1, aes(sample = fitted))+
				stat_qq(size=1,pch=21)+
				labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
				stat_qq_line(linetype=2, col='red')+
				theme_bw()

			diagnostics_fit1 <- resVfit_fit1+res_fit1_hist+res_fit1_qq
			
			ggsave(diagnostics_fit1, file = 'EDNA/data_edna/brms/figures/diagnostics_model1_gaussian_log10data.png', device = 'png', unit = 'in', height = 4, width = 8, dpi = 850)	





























