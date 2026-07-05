# plot tree mortality
# SCE 12 Jan 2023
#  Credit for basic workflow for downloading/top subsetting/cleaning from
# C. Bueno de mesquita script 

# -- SETUP -----
# clean up enviro, read in needed libraries
rm(list=ls())
#library(readxl)
library(tidyverse)
#library(vegan)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, NaN, "NA", "NaN", ".")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utility_functions/utility_functions_all.R")


#cliff function summarySE

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#TODO 


# download data -----------------------------------------------------------


# Tree Mortality by PC1
# 1982 - 2016
# by Cliff Bueno de Mesquita
# -- SETUP -----
# clean up enviro, read in needed libraries
# rm(list=ls())
# library(RCurl)
# script <- getURL("https://raw.githubusercontent.com/cliffbueno/Functions/master/Summary.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
# library(readxl)
# library(tidyverse)
# library(vegan)
# library(naniar)
# library(FSA)
# library(plyr)
# options(stringsAsFactors = F)
# theme_set(theme_bw())
# na_vals <- c(" ", "", NA, NaN, "NA", "NaN", ".")

# -- FUNCTIONS -----
# also just run the whole utilities R code
# set up functions to read in tabular datasets from EDI dynamically
# SCE + CTW code to determine most recent version of package ID and read in current dataset on EDI
#function to determine current version of data package on EDI
## ----- UTILITY FUNCTIONS ----
## Need to improve documentation and consider splitting into multiple files?
findNonNumeric<-function(x){
  unique(suppressWarnings(x[is.na(as.numeric(x))]))
}

# #borrowed from metajam
# tabularize_eml <- function(eml, full = FALSE) {
#   
#   if (any(class(eml) == "emld")) {
#     eml <- eml
#   } else if (is.character(eml) | is.raw(eml)) {
#     eml <- emld::as_emld(eml)
#   } else {
#     stop("The EML input could not be parsed.")
#   }
#   
#   metadata <- eml %>%
#     unlist() %>%
#     tibble::enframe()
#   
#   if (full == FALSE) {
#     metadata <- metadata %>%
#       dplyr::mutate(name = dplyr::case_when(
#         grepl("schemaLocation", name) ~ "eml.version",
#         grepl("title", name) ~ "title",
#         grepl("individualName", name) ~ "people",
#         grepl("abstract", name) ~ "abstract",
#         grepl("keyword", name) ~ "keyword",
#         grepl("geographicDescription", name) ~ "geographicCoverage.geographicDescription",
#         grepl("westBoundingCoordinate", name) ~ "geographicCoverage.westBoundingCoordinate",
#         grepl("eastBoundingCoordinate", name) ~ "geographicCoverage.eastBoundingCoordinate",
#         grepl("northBoundingCoordinate", name) ~ "geographicCoverage.northBoundingCoordinate",
#         grepl("southBoundingCoordinate", name) ~ "geographicCoverage.southBoundingCoordinate",
#         grepl("beginDate", name) ~ "temporalCoverage.beginDate",
#         grepl("endDate", name) ~ "temporalCoverage.endDate",
#         grepl("taxonRankValue", name) ~ "taxonomicCoverage",
#         grepl("methods", name) ~ "methods",
#         grepl("objectName", name) ~ "objectName",
#         grepl("online.url", name) ~ "url"
#       )) %>%
#       dplyr::filter(!is.na(name)) %>%
#       dplyr::mutate(value = stringr::str_trim(value)) %>%
#       dplyr::distinct() %>%
#       dplyr::group_by(name) %>%
#       dplyr::summarize(value = paste(value, collapse = "; ")) %>%
#       dplyr::mutate(value = gsub("\n", "", value)) #without this, fields get truncated in Excel
#   }
#   
#   return(metadata)
# }
# 
# #function to determine current version of data package on EDI
# getCurrentVersion<-function(edi_id){
#   require(magrittr)
#   versions=readLines(paste0('https://pasta.lternet.edu/package/eml/knb-lter-nwt/', edi_id), warn=FALSE)%>%
#     as.numeric()%>%(max)
#   packageid=paste0('knb-lter-nwt.', edi_id, '.', versions)
#   return (packageid)
# }
# 
# #function to download the EML file from EDI
# getEML<-function(packageid){
#   require(magrittr)
#   myurl<-paste0("https://portal.lternet.edu/nis/metadataviewer?packageid=",
#                 packageid,
#                 "&contentType=application/xml")
#   #myeml<-xml2::download_html(myurl)%>%xml2::read_xml()%>%EML::read_eml()
#   myeml<-xml2::read_xml(paste0("https://portal.lternet.edu/nis/metadataviewer?packageid=",
#                                packageid,
#                                "&contentType=application/xml"))%>%EML::read_eml()
# }
# 
# #function to get a single element anywhere in an eml
# eml_get_simple <- function(x, element, from = "list", ...){
#   doc <- as.character(emld::as_json(emld::as_emld(x, from = from)))
#   out <- jqr::jq(doc, paste0("..|.", element, "? // empty"))
#   json <- jqr::combine(out)
#   robj <- jsonlite::fromJSON(json, simplifyVector = FALSE)
#   return(robj)
# }
# 
# #function for plot aesthetics
# use_theme = function(){ 
#   theme_bw()+
#     theme(
#       #this is size of the font on the panels
#       strip.text.y = element_text(size = 3),
#       #yaxis font size
#       axis.text.y = element_text(size = 4),
#       #axis.title.y = element_blank(),
#       axis.text.x=element_text(angle=-90)
#     )
# }
# 
# # function to determine data package version number only (not return full data package ID)
# getPackageVersion<-function(edi_id, site = "nwt"){
#   versions=readLines(paste0('https://pasta.lternet.edu/package/eml/knb-lter-', site, '/', edi_id), warn=FALSE)
#   currentV <- max(as.numeric(versions))
#   return(currentV)
# }
# 
# # function to get entity ID for current data package version
# getEntityId <- function(edi_id, version, site = "nwt", datanum = 1){
#   entID <- readLines(paste0('https://pasta.lternet.edu/package/eml/knb-lter-', site, '/', edi_id, "/", version, "/"), warn=FALSE)[datanum]
#   entID <- gsub(paste0("http.*/",edi_id,"/",version,"/"), "", entID) # remove all chars except what comes after last /
#   return(entID)
# }
# 
# # function to read in tabular csv dataset for data package (could make more generic with read table, but should know what you're reading in to use)
# getTabular <- function(edi_id, na_vals = c("", "NA", NA, NaN, ".", "NaN", " "), site = "nwt", datanum = 1){
#   require(readr)
#   v <- getPackageVersion(edi_id, site = site)
#   id <- getEntityId(edi_id, v, site = site, datanum = datanum)
#   dat <- readr::read_csv(paste0("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-", site, ".", edi_id, ".", v, 
#                                 "&entityid=", id),
#                          trim_ws =TRUE, na = na_vals)
#   dat <- as.data.frame(dat)
#   print(paste0("Reading in knb-lter-", site, ".", edi_id, ".", v))
#   return(dat)
# }

# -- GET DATA -----
# Tom Veblen Tree Dataset
PP_plot_data <- getTabular(207, datanum = 1)
PP_tree_data <- getTabular(207, datanum = 2)
PP_seedling_data <- getTabular(207, datanum = 3)


# let's just look at the raw and try to remember how this is structured

View (PP_tree_data %>%
        filter(plot == 'BW3'))

# appears you can only tell for sure what was censused
# when if there's something that died in that period...

#adding in topographic position from
#https://www.sciencedirect.com/science/article/pii/S0378112714007476#t0005

extra_plot_info=data.frame(
  plot=c('BW2', 'BW3','MRS4', 'MRS5', 'BL6', 'MRS7', 'MRS4_gap', 'BL6_gap', 'MRS1','MRS8-9-10', 'MRS1_gap'),
  topo=c('Xeric', 'Mesic', 'Xeric', 'Hydric', 'Mesic', 'Xeric', 'Xeric/Mesic' , 'Mesic', 'Xeric', 'Xeric/Mesic', 'Xeric'),
  age= c('old', 'old', 'old', 'old', 'old', 'old', 'old', 'young', 'young', 'young', 'young')
)

# Gap plots 6-10 and 12-21 are over a broad area west of near MRS4 (old, xeric/mesic)
# Gap plots 11 and 22-30 are near BL6 (old, xeric/mesic) (old, mesic)
# MRS8, 9 & 10 (xeric/mesic, young)
# Gap plots 1 through 5 are in a tight cluster on the north side of the ski trail north of MRS1 (young, xeric)

PP_tree_data=PP_tree_data%>%dplyr::mutate(
  plot_group=plot)%>%
  dplyr::mutate(
    plot_group=replace(plot_group, plot_group%in%c('Gap1',  'Gap2',  'Gap3',  'Gap4',  'Gap5'), 'MRS1_gap'))%>%
  dplyr::mutate(
    plot_group=replace(plot_group, plot_group%in%c('MRS8',  'MRS9',  'MRS10'), 'MRS8-9-10'))%>%
  dplyr::mutate(plot_group=replace(plot_group, plot_group%in%c('Gap6',  'Gap7',  'Gap8',  'Gap9',  'Gap10',
                                                               'Gap12',  'Gap13',  'Gap14',  'Gap15',  'Gap15',
                                                               'Gap16',  'Gap17',  'Gap18',  'Gap19',  'Gap20', 'Gap21'), 'MRS4_gap'))%>%
  dplyr::mutate(plot_group=replace(plot_group, plot_group%in%c('Gap11',  'Gap22',  'Gap23',  'Gap24',  'Gap25',
                                                               'Gap26',  'Gap27',  'Gap28',  'Gap29',  'Gap30'), 'BL6_gap'))%>%
  dplyr::left_join(., PP_plot_data %>%dplyr::select(plot, install_year))

#fill in more missing MD from paper
PP_tree_data$install_year[PP_tree_data$plot_group=='MRS1_gap']=1983
PP_tree_data$install_year[PP_tree_data$plot_group=='BL6_gap']=1983
PP_tree_data$install_year[PP_tree_data$plot_group=='MRS4_gap']=1983

#start of dead period for things with no dp_start assumed to be the install year
#probably things with no start AND end to the dp were dead at the initial survey??
#
PP_tree_data$dp_start[is.na(PP_tree_data$dp_start)&PP_tree_data$dead==1&!is.na(PP_tree_data$dp_end)]=
  PP_tree_data$install_year[is.na(PP_tree_data$dp_start)&PP_tree_data$dead==1&!is.na(PP_tree_data$dp_end)]

PP_tree_data$duration=(PP_tree_data$dp_end-PP_tree_data$dp_start)+1

##First plot height class by measurement year
num_vars<-PP_tree_data %>%
  dplyr::select_if(., is.numeric)  %>%
  names()

PP_tree_data$dbhmean <- rowMeans(PP_tree_data[c('dbh1', 'dbh3')], na.rm=TRUE)

ggplot(PP_tree_data, aes(x=yrcht, y=dbhmean, group=`spec`))+ geom_point(aes(col=spec))

##Plot percent mortality for each death period and by height class
PP_tree_data = PP_tree_data%>%tidyr::unite(.,  period, dp_start, dp_end, sep="-" )

#confirmed, one alive at end fo period, 7 died during period
# so Cliff's mortality calcs need adjusting
# so that mortality is dead/(dead+alive)
View (PP_tree_data %>%
        filter(plot == 'BW2'&spec == 'POTR'))
#unique (PP_tree_data[,c('plot', 'period')])

#note in the published data, they calculated
# percent mortality based on the number of trees at plot
# installation but I think properly it should be done off
# those at the beginning of the period

#assign a status for all trees at each period
# dp = died previous period
# d = died this period
# a = alive this period
# There was a problem with the tree. column, rename


tst<-tidyr::expand(PP_tree_data, tidyr::nesting(plot, tree.), period)

#so far assumes there was a BW3 census in this period
View (tst %>%filter(plot == 'BW3' & period == '2014-2016'))

#lots of periods here, seems it has been censussed tons of times!
unique (PP_tree_data %>% filter(plot == 'BW3') %>%dplyr::pull(period))


names(PP_tree_data)[2] <- c("tree.")
mortality=tidyr::expand(PP_tree_data, tidyr::nesting(plot, tree.), period)%>%
  dplyr::left_join(., PP_tree_data%>%dplyr::select(plot, tree., spec, period,
                                                   dead, hc1,biotic_attack, biotic_agent, plot_group)%>%
                     dplyr::rename(mort_period=period))%>%
  dplyr::arrange(plot, tree., period)

#this is somewhat weird because there is both 1993-2007 and also 1995-2007
# how is this possible?
# does not exist in the raw data but instead is an error in the nesting 
# statement above.

mortnew=tidyr::expand(PP_tree_data, tidyr::nesting(plot, tree.), period) 
nrow (mortnew)

#get rid of censuses that don't exist (SCE thinks this is right 
# but we would need to double check w Robbie!)
mortnew = mortality %>%
  inner_join(., PP_tree_data %>%
               select(plot, period) %>%
               distinct())

#this then gets rid of a ton of faux tree censuses
nrow (mortnew)

mortality <- mortnew %>%
  dplyr::left_join(., PP_tree_data%>%dplyr::select(plot, tree., spec, period,
                                                   dead, hc1,biotic_attack, biotic_agent, plot_group)%>%
                     dplyr::rename(mort_period=period))%>%
  dplyr::arrange(plot, tree., period)
  


# unique (mortality %>% filter(plot == 'BW3') %>%dplyr::pull(period))
# 
# View (PP_tree_data %>%
#         select(plot, period) %>%
#         distinct() %>%
#         arrange(plot, period))



#but, ignoring that for now because the very long intercensus period
# is dropped anyhow...

mortality$status=NA
for (i in 1:nrow(mortality)){
  if (mortality$mort_period[i]=='NA-NA'){
    mortality$status[i]='a'
  }else if (mortality$mort_period[i]==mortality$period[i]){
    mortality$status[i]='d'
  }else if (mortality$status[i-1]%in%(c('d', 'dp'))&mortality$tree.[i]==mortality$tree.[i-1]){
    mortality$status[i]='dp'
  }else{
    mortality$status[i]='a'
  }
}

#SCE let's back out for the period 2014-2016 and PIEN and BW3, what happened
# to give a d = NA 

#seems pretty clear that nothing died in this period
View (mortality %>%
        filter(spec == 'PIEN'& plot_group == 'BW3'))
      
#but not looking better after fixing the mortnew above
# there is not a fake period in here, continuing w the code...


#summarize percent mortality by plot and period over all species
mortrates=mortality%>%
  dplyr::group_by(plot_group, period, status)%>%
  dplyr::summarise(n=dplyr::n())%>%
  tidyr::spread(., status, n)%>%
  dplyr::filter(., !period%in%c('1982-1982', '1983-1983', '1986-1986', '2016-2016', 'NA-NA'))%>%
  dplyr::left_join(., PP_plot_data, by=c('plot_group'='plot'))

#more from the ms - install years for the groups
mortrates$install_year[mortrates$plot_group=="BL6_gap"]=1982
mortrates$install_year[mortrates$plot_group=="MRS4_gap"]=1983
mortrates$install_year[mortrates$plot_group=="MRS8-9-10"]=1986
mortrates$install_year[mortrates$plot_group=="MRS1_gap"]=1983

#remove census period prior to install year
mortrates=mortrates %>%
  dplyr::filter(!(install_year==1986&period=='1984-1986'))%>%
  dplyr::left_join(., unique(PP_tree_data[, c('period', 'duration')]))%>%
  dplyr::left_join(., extra_plot_info, by=c('plot_group'='plot'))%>%
  dplyr::mutate(perc_mort=d/a)%>%
  dplyr::mutate(perc_mort_ann=perc_mort/duration)

#calculate mean ext summer over the periods
# clims=mortrates%>%dplyr::ungroup()%>%
#   dplyr::select(period)%>%
#   distinct()%>%
#   tidyr::separate(., period, into=c('dp_start', 'dp_end'))
# 
# clims$mean_ext_summer=NA
# for (j in 1:nrow (clims)){
#   clims$mean_ext_summer[j]=mean(extsummer$sumallPC1[extsummer$eco_year>=clims$dp_start[j]&
#                                                       extsummer$eco_year<=clims$dp_end[j]])
# }
# 
# mortrates=mortrates%>%dplyr::left_join(.,
#                                        clims%>%tidyr::unite(.,  period, dp_start, dp_end, sep="-" ))


topocolor=mortrates%>%
  ungroup()%>%
  dplyr::select(plot_group, topo)%>%
  distinct()%>%
  dplyr::mutate(
    color=case_when(
      topo=='Xeric' ~ 'brown',
      topo=='Xeric/Mesic' ~ 'olivedrab4',
      topo=='Mesic' ~ 'green',
      topo=='Hydric' ~ 'blue'))%>%
  dplyr::arrange(topo)

topocolor=dplyr::bind_rows(
  topocolor%>%dplyr::filter(topo=='Hydric'),
  topocolor%>%dplyr::filter(topo=='Mesic'),
  topocolor%>%dplyr::filter(topo=='Xeric/Mesic'),
  topocolor%>%dplyr::filter(topo=='Xeric'))

#reorder by moisture
mortrates$plot_group <- factor(mortrates$plot_group, levels=topocolor$plot_group)

# Calculate actual percent
mortrates$PercAnnMort <- mortrates$perc_mort_ann*100
# Can't use 15 yr periods
mortrates_sub <- subset(mortrates, duration < 4)
# Make factor to summarize by
# mortrates_sub$extsumfact <- as.factor(mortrates_sub$mean_ext_summer)
# meanmort <- summarySE(mortrates_sub, measurevar ="PercAnnMort",groupvars=c("extsumfact","topo"))
# meanmort$char <- as.character(meanmort$extsumfact)
# meanmort$PC1 <- as.numeric(meanmort$char)


meanmort <- summarySE(mortrates_sub, measurevar ="PercAnnMort",groupvars=c("period","topo"))


# Plot by moisture
ggplot(meanmort, aes(x=period, y=PercAnnMort, colour = topo)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin=PercAnnMort-se,ymax=PercAnnMort+se),width=0,alpha = 0.5) +
  geom_smooth(method = lm, se = F) +
  labs(x = "period",
       y = "% Annual Tree Mortality",
       colour = "Site Type") +
    theme(legend.position = c(0.14, 0.85),
          legend.background = element_rect(color = "black"),
          legend.title = element_text(face="bold", size = 14),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1, "line"), 
          axis.text = element_text(size = 16), 
          axis.title = element_text(face="bold",size=18))

# # Plot all
# ggplot(meanmort, aes(x=PC1, y=PercAnnMort)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin=PercAnnMort-se,ymax=PercAnnMort+se),width=0) +
#   geom_smooth(method = lm, se = F) +
#   labs(x = "PC1",
#        y = "% Annual Tree Mortality") +
#   theme(axis.text = element_text(size = 16), 
#         axis.title = element_text(face="bold",size=18))
# 
# # Regressions on Means
# hydm <- lm(PercAnnMort ~ PC1, data = subset(meanmort, topo == "Hydric"))
# summary(hydm) # NS
# mesm <- lm(PercAnnMort ~ PC1, data = subset(meanmort, topo == "Mesic"))
# summary(mesm) # NS
# xerm <- lm(PercAnnMort ~ PC1, data = subset(meanmort, topo == "Xeric"))
# summary(xerm) # Only two points
# xemm <- lm(PercAnnMort ~ PC1, data = subset(meanmort, topo == "Xeric/Mesic"))
# summary(xemm) # NS
# allm <- lm(PercAnnMort ~ PC1, data = meanmort)
# summary(allm) # R2 = 0.24, p = 0.056
# 
# # Regressions on Whole Data
# gm <- lm(PercAnnMort ~ topo+mean_ext_summer, data = mortrates_sub)
# summary(gm) # Topo and PC1 significant
# hydm1 <- lm(PercAnnMort ~ mean_ext_summer, data = subset(mortrates_sub, topo == "Hydric"))
# summary(hydm1) # NS
# mesm1 <- lm(PercAnnMort ~ mean_ext_summer, data = subset(mortrates_sub, topo == "Mesic"))
# summary(mesm1) # R2 = 0.22, p = 0.056
# xerm1 <- lm(PercAnnMort ~ mean_ext_summer, data = subset(mortrates_sub, topo == "Xeric"))
# summary(xerm1) # R2 = 0.31, p = 0.002
# xemm1 <- lm(PercAnnMort ~ mean_ext_summer, data = subset(mortrates_sub, topo == "Xeric/Mesic"))
# summary(xemm1) # NS
# allm1 <- lm(PercAnnMort ~ mean_ext_summer, data = mortrates_sub)
# summary(allm1) # R2 = 0.15, p = 0.002
# 

################################## By Species ##############################################
#summarize percent mortality by plot and period
mortrates=mortality%>%
  dplyr::group_by(plot_group, period, status, spec)%>%
  dplyr::summarise(n=dplyr::n())%>%
  tidyr::spread(., status, n)%>%
  dplyr::filter(., !period%in%c('1982-1982', '1983-1983', '1986-1986','2016-2016', 'NA-NA'))%>%
  dplyr::left_join(., PP_plot_data, by=c('plot_group'='plot'))



#SCE added this line I am pretty sure if there are a but not d then d is 
# actually 0 - missing from Cliff's logic so need to double check

View (mortrates %>%
        filter(is.na(d)))

#sce infilled these w 0
tst = mortrates %>%
  mutate(d_ch = ifelse(is.na(d)&!is.na(a), 0, d))

#View (tst %>% select(d_ch, everything())) # seems right?

#so change the d calc to fill in the 0s
mortrates <- mortrates %>%
  mutate(d = ifelse(is.na(d)&!is.na(a), 0, d))


#more from the ms - install years for the groups
mortrates$install_year[mortrates$plot_group=="BL6_gap"]=1982
mortrates$install_year[mortrates$plot_group=="MRS4_gap"]=1983
mortrates$install_year[mortrates$plot_group=="MRS8-9-10"]=1986
mortrates$install_year[mortrates$plot_group=="MRS1_gap"]=1983

#remove census period prior to install year
mortrates=mortrates %>%
  dplyr::filter(!(install_year==1986&period=='1984-1986'))%>%
  dplyr::left_join(., unique(PP_tree_data[, c('period', 'duration')]))%>%
  dplyr::left_join(., extra_plot_info, by=c('plot_group'='plot'))%>%
  #sce changed this line I think it is not right
  #dplyr::mutate(perc_mort=d/a)%>%
  dplyr::mutate(perc_mort=d/(d+a))%>%
  dplyr::mutate(perc_mort_ann=perc_mort/duration)

#mortrates=mortrates%>%dplyr::left_join(.,
#                                       clims%>%tidyr::unite(.,  period, dp_start, dp_end, sep="-" ))


#reorder by moisture
mortrates$plot_group <- factor(mortrates$plot_group, levels=topocolor$plot_group)

# Can't use 15 yr periods
#mortrates$extsumfact <- as.factor(mortrates$mean_ext_summer)
mortrates$PercAnnMort <- mortrates$perc_mort_ann*100
mortrates_sub <- subset(mortrates, duration < 4)

#sanity test this now looks good
# mortality rates all less than 0
# and are also NA if nothing alive
View (mortrates_sub %>%
        filter(spec == 'POTR'))

#as we have it we have removed the 100%
# mortality from a few plots where the single last tree died
# not sure if this makes sense or not
View (mortrates %>%
        filter(is.na(PercAnnMort)))

mortrates_all = mortrates %>%
  mutate(PercAnnMort_adj = ifelse((is.na(a)&!is.na(d)),
         1, PercAnnMort))

# View (mortrates_all %>%
#         filter(is.na(PercAnnMort)))

# Summarize
#meanmortsp <- summarySE(mortrates_sub, measurevar ="PercAnnMort",groupvars=c("extsumfact","topo","spec"))
meanmortsp <- summarySE(mortrates_sub, measurevar ="PercAnnMort",groupvars=c("period","topo","spec"))
#meanmortsp$char <- as.character(meanmortsp$extsumfact)
#meanmortsp$PC1 <- as.numeric(meanmortsp$char)
meanmortsp <- subset(meanmortsp, PercAnnMort != "NA")
meanmortsp$spec <- as.factor(meanmortsp$spec)
levels(meanmortsp$spec)
table(meanmortsp$spec) # Remove PIFL
meanmortsp <- subset(meanmortsp, spec != "PIFL")

# Plot with moisture color
facet_names <-  c(`ABLA` = "Subalpine Fir",`PICO` = "Lodgepole Pine",`PIEN` = "Engelmann Spruce",`POTR` = "Quaking Aspen")


#something definitely not quite right
# as there is an aspen plot with >100% mortality?
# but then again maybe it wa
ggplot(meanmortsp, aes(x=period, y=PercAnnMort, colour = topo)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin=PercAnnMort-se,ymax=PercAnnMort+se),width=0,alpha = 0.5) +
  #geom_smooth(method = lm, se = F) +
  labs(x = "period",
       y = "% Annual Tree Mortality",
       colour = "Site Type") +
  facet_wrap(~ spec, scales = "free_y", labeller = as_labeller(facet_names)) +
  theme(legend.position = "right",
        legend.background = element_rect(color = "black"),
        legend.title = element_text(face="bold", size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "line"), 
        axis.text = element_text(size = 16), 
        axis.title = element_text(face="bold",size=18))

# 
# View (meanmortsp %>%
#         filter(spec == 'POTR'))

# # Plot
# facet_names <-  c(`ABLA` = "Subalpine Fir",`PICO` = "Lodgepole Pine",`PIEN` = "Engelmann Spruce",`POTR` = "Quaking Aspen")
# ggplot(meanmortsp, aes(x=PC1, y=PercAnnMort, group = spec)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin=PercAnnMort-se,ymax=PercAnnMort+se),width=0) +
#   geom_smooth(method = lm, se = F) +
#   labs(x = "PC1",
#        y = "% Annual Tree Mortality") +
#   facet_wrap(~ spec, scales = "free_y", labeller = as_labeller(facet_names)) +
#   theme(axis.text = element_text(size = 16), 
#         axis.title = element_text(face="bold",size=18),
#         strip.text = element_text(size = 14))
# 
# # Regressions on Means
# ablam <- lm(PercAnnMort ~ PC1, data = subset(meanmortsp, spec == "ABLA"))
# summary(ablam) # NS
# picom <- lm(PercAnnMort ~ PC1, data = subset(meanmortsp, spec == "PICO"))
# summary(picom) # NS
# pienm <- lm(PercAnnMort ~ PC1, data = subset(meanmortsp, spec == "PIEN"))
# summary(pienm) # NS
# potrm <- lm(PercAnnMort ~ PC1, data = subset(meanmortsp, spec == "POTR"))
# summary(potrm) # R2 = 0.99, p = 0.048
# 
# # Regressions on whole data
# gm1 <- lm(perc_mort_ann ~ spec+topo+mean_ext_summer, data = mortrates_sub)
# summary(gm1) # Species significant, moisture and PC1 not.
# ablam1 <- lm(perc_mort_ann ~ mean_ext_summer, data = subset(mortrates_sub, spec == "ABLA"))
# summary(ablam1) # R2 = 0.15, p = 0.02
# picom1 <- lm(perc_mort_ann ~ mean_ext_summer, data = subset(mortrates_sub, spec == "PICO"))
# summary(picom1) # NS
# pienm1 <- lm(perc_mort_ann ~ mean_ext_summer, data = subset(mortrates_sub, spec == "PIEN"))
# summary(pienm1) # NS
# potrm1 <- lm(perc_mort_ann ~ mean_ext_summer, data = subset(mortrates_sub, spec == "POTR"))
# summary(potrm1) # R2 = 0.99, p = 0.048
# 
# # Most accurate graph with ABLA and POTR increasing with PC1
# ggplot(meanmortsp, aes(x=PC1, y=PercAnnMort, colour = spec)) +
#   geom_point(size = 3, alpha = 0.5) +
#   geom_errorbar(aes(ymin=PercAnnMort-se,ymax=PercAnnMort+se),width=0,alpha = 0.5) +
#   labs(x = "Extended Summer (PC1 Axis Score)",
#        y = "% Annual Tree Mortality",
#        colour = NULL) +
#   geom_smooth(method = lm, se = F, data = subset(meanmortsp, spec == "ABLA" | spec == "POTR")) +
#   theme(legend.position = c(0.15,0.75),
#         legend.background = element_rect(color = "black"),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 12),
#         legend.key.size = unit(2, "line"), 
#         axis.text = element_text(size = 16), 
#         axis.title = element_text(face="bold",size=18))

# Cliff's notes
# Update after talking to Robbie
# Too few POTR trees to have confidence so remove
# PICO mortality largely due to density dependence
# Just look at ABLA and PIEN and in plots were Robbie says climate could play a role instead of density dependence
# ReSummarize
# Summarize
#meanmortsp <- summarySE(mortrates_sub, measurevar ="PercAnnMort",groupvars=c("extsumfact","topo","spec"))
#meanmortsp <- summarySE(mortrates_sub, measurevar ="PercAnnMort",groupvars=c("period","topo","spec"))
#meanmortsp$char <- as.character(meanmortsp$extsumfact)
#meanmortsp$PC1 <- as.numeric(meanmortsp$char)
# meanmortsp <- subset(meanmortsp, PercAnnMort != "NA")
# meanmortsp$spec <- as.factor(meanmortsp$spec)
# levels(meanmortsp$spec)
# table(meanmortsp$spec) # Remove PIFL
# meanmortsp <- subset(meanmortsp, spec != "PIFL")


# Graph with Robbie recommended Species and Plots
# SCE note I am unsure why Robbie suggested only these plots
# I think there are too few aspen and the spruce apparently have density 
# dependence
# just following Cliff's notes here
mortrates_sub <- subset(mortrates_sub,plot_group == "BW3"|plot_group=="MRS4"|plot_group=="MRS5"|plot_group=="BL6"|plot_group=="MRS7")

#add back in the 100% mortality
mortrates_sub_all = mortrates_sub  %>%
  mutate(PercAnnMort_adj = ifelse((is.na(a)&!is.na(d)),
                                  1, PercAnnMort))

# meanmortsp <- summarySE(mortrates_sub, measurevar ="PercAnnMort",groupvars=c("period","plot_group","spec"))
# #meanmortsp$char <- as.character(meanmortsp$extsumfact)
# #meanmortsp$PC1 <- as.numeric(meanmortsp$char)
# meanmortsp <- subset(meanmortsp, PercAnnMort != "NA")
# meanmortsp$spec <- as.factor(meanmortsp$spec)
# levels(meanmortsp$spec)
# table(meanmortsp$spec) # Remove PIFL
# meanmortsp <- subset(meanmortsp, spec == "ABLA" | spec == "PIEN")
# meanmortsp <- subset(meanmortsp, plot_group == "BW3"|plot_group=="MRS4"|plot_group=="MRS5"|plot_group=="BL6"|plot_group=="MRS7")
# mortrates_sub <- subset(mortrates_sub,plot_group == "BW3"|plot_group=="MRS4"|plot_group=="MRS5"|plot_group=="BL6"|plot_group=="MRS7")
# Regressions
# ablam1 <- lm(perc_mort_ann ~ mean_ext_summer, data = subset(mortrates_sub, spec == "ABLA"))
# summary(ablam1) # R2 = 0.23, p = 0.02
# pienm1 <- lm(perc_mort_ann ~ mean_ext_summer, data = subset(mortrates_sub, spec == "PIEN"))
# summary(pienm1) # NS
# Graph with Robbie recommended Species and Plots
# Significant line for ABLA

#sce compared this to cliff's figure and seems about like he had intended
# eg one point per plot per period (with Cliff's axis being extended summer)
# as that was the framework then
ggplot(meanmortsp %>%
         filter(spec %in% c('PIEN', 'ABLA')), aes(x=period, y=PercAnnMort, colour = spec)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin=PercAnnMort-se,ymax=PercAnnMort+se),width=0,alpha = 0.5) +
  #labs(x = "Extended Summer (PC1 Axis Score)",
  labs(x = "Census period",
       y = "% Annual Tree Mortality",
       colour = NULL) +
  scale_colour_discrete(labels = c("Subalpine Fir","Engelmann Spruce")) +
  geom_smooth(method = lm, se = F, data = subset(meanmortsp, spec == "ABLA")) +
  theme(legend.position = c(0.20,0.75),
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(2, "line"), 
        axis.text = element_text(size = 16), 
        axis.title = element_text(face="bold",size=18))


#sce remaking by hydric/xeric and period
meanmortsp <- summarySE(mortrates_sub %>%
                          filter(plot_group == "BW3"|plot_group=="MRS4"|plot_group=="MRS5"|plot_group=="BL6"|plot_group=="MRS7"), measurevar ="PercAnnMort",groupvars=c("period","topo","spec"))
#meanmortsp$char <- as.character(meanmortsp$extsumfact)
#meanmortsp$PC1 <- as.numeric(meanmortsp$char)
meanmortsp <- subset(meanmortsp, PercAnnMort != "NA")
meanmortsp$spec <- as.factor(meanmortsp$spec)
levels(meanmortsp$spec)
table(meanmortsp$spec) # Remove PIFL
meanmortsp <- subset(meanmortsp, spec == "ABLA" | spec == "PIEN")
#meanmortsp <- subset(meanmortsp, plot_group == "BW3"|plot_group=="MRS4"|plot_group=="MRS5"|plot_group=="BL6"|plot_group=="MRS7")
#mortrates_sub <- subset(mortrates_sub,plot_group == "BW3"|plot_group=="MRS4"|plot_group=="MRS5"|plot_group=="BL6"|plot_group=="MRS7")
# Regressions
# ablam1 <- lm(perc_mort_ann ~ mean_ext_summer, data = subset(mortrates_sub, spec == "ABLA"))
# summary(ablam1) # R2 = 0.23, p = 0.02
# pienm1 <- lm(perc_mort_ann ~ mean_ext_summer, data = subset(mortrates_sub, spec == "PIEN"))
# summary(pienm1) # NS
# Graph with Robbie recommended Species and Plots
# Significant line for ABLA



#sce compared this to cliff's figure and seems about like he had intended
# eg one point per plot per period (with Cliff's axis being extended summer)
# as that was the framework then
g1<- ggplot(meanmortsp, aes(x=period, y=PercAnnMort, colour = spec)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin=PercAnnMort-se,ymax=PercAnnMort+se),width=0,alpha = 0.5) +
  #labs(x = "Extended Summer (PC1 Axis Score)",
  labs(x = "Census period",
       y = "% Annual Tree Mortality",
       colour = NULL) +
  scale_colour_discrete(labels = c("Subalpine Fir","Engelmann Spruce")) +
  geom_smooth(method = lm, se = F, data = subset(meanmortsp, spec == "ABLA")) +
  theme(legend.position = c(0.20,0.75),
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(2, "line"), 
        axis.text = element_text(size = 16), 
        axis.title = element_text(face="bold",size=18))+
  facet_wrap(~topo)+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

ggsave(g1, file = 'plots/mortality_by_census.jpg', height =4, width =12)


#other option, just make a weighted mean over all plots per period by top
meanmortsp_adj =
  mortrates_sub_all %>%
  filter(plot_group == "BW3"|plot_group=="MRS4"|plot_group=="MRS5"|plot_group=="BL6"|plot_group=="MRS7")%>%
  mutate(ct_a = ifelse(!is.na(a), a, 0),
         ct_d = ifelse(!is.na(d), d, 0),
         sum_tree = ct_a+ct_d) %>%
  ungroup() %>%
  filter(sum_tree!=0)%>%
  dplyr::group_by(topo, spec,period) %>%
  #dplyr::summarize(PercAnnMort = mean(PercAnnMort_adj))
  dplyr::summarize(PercAnnMort = weighted.mean(PercAnnMort_adj, sum_tree), .groups = 'drop') %>%
  tidyr::separate(., col = 'period', into = c('start', 'end'), sep = '-') %>%
  mutate(start = as.numeric(start), end = as.numeric(end),
         year = (start+end)/2) %>%
  filter(spec %in% c('ABLA', 'PIEN')) %>%
  dplyr::group_by(spec, topo) %>%
  dplyr::mutate(mean_mort = mean (PercAnnMort,na.rm = TRUE), .groups = 'drop') %>%
  mutate(pos_neg = ifelse(PercAnnMort> mean_mort, 'pos', 'neg'))



g2<- ggplot(meanmortsp, aes(x=period, y=PercAnnMort, colour = spec)) +
  geom_col()+
  #geom_point(size = 3, alpha = 0.5) +
  #geom_errorbar(aes(ymin=PercAnnMort-se,ymax=PercAnnMort+se),width=0,alpha = 0.5) +
  #labs(x = "Extended Summer (PC1 Axis Score)",
  labs(x = "Census period",
       y = "% Annual Tree Mortality",
       colour = NULL) +
  #scale_colour_discrete(labels = c("Subalpine Fir","Engelmann Spruce")) +
  theme(legend.position = c(0.20,0.75),
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(2, "line"), 
        axis.text = element_text(size = 16), 
        axis.title = element_text(face="bold",size=18))+
  facet_grid(spec~topo)+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))



# Plot with moisture color
facet_names <-  c(`ABLA` = "Subalpine Fir",
                  `PICO` = "Lodgepole Pine",
                  `PIEN` = "Engelmann Spruce",
                  `POTR` = "Quaking Aspen",
                  `Hydric` = "Hydric",
                  `Mesic` = "Mesic",
                  `Xeric` = "Xeric")

g3<- ggplot(meanmortsp_adj %>%
              filter(spec %in% c('ABLA', 'PIEN')), aes(x=year, y=PercAnnMort, fill = pos_neg)) +
  geom_col()+
  #geom_point(size = 3, alpha = 0.5) +
  #geom_errorbar(aes(ymin=PercAnnMort-se,ymax=PercAnnMort+se),width=0,alpha = 0.5) +
  #labs(x = "Extended Summer (PC1 Axis Score)",
  labs(x = "Census period",
       y = "% Annual Tree Mortality",
       colour = NULL) +
  scale_fill_manual(values = c("green4", "chocolate4")) +
  #scale_colour_discrete(labels = c("Subalpine Fir","Engelmann Spruce")) +
  # theme(legend.position = c(0.20,0.75),
  #       legend.background = element_rect(color = "black"),
  #       legend.title = element_blank(),
  #       legend.text = element_text(size = 12),
  #       legend.key.size = unit(2, "line"), 
  #       axis.text = element_text(size = 16), 
  #       axis.title = element_text(face="bold",size=18))+
  facet_grid(spec~topo, labeller = as_labeller(facet_names))+
  theme(legend.position = "none")
  


ggsave(g3, file = 'plots/mortality_by_census_by_spp.jpg', height =8, width =10,
       scale = 0.5)




  



measurevar ="PercAnnMort",groupvars=c("period","topo","spec"))


