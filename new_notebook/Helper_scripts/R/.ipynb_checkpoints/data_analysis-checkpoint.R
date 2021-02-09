### DATA_ANALYSIS.R

#define colour palette to be used
cbPalette <- c("#009E73", "#F0E442", "#0072B2")

###
#Read in data
###

#read in OTU Resp
OTU_R <- read_csv("../data/Francisca/respirationandODData//RespDataLinear_OTUs_evNev.csv",col_types = cols())

OTU_R <- OTU_R %>% mutate(OTU = as.numeric(str_replace(Sps,"R_","")),
       Resp = -Resp,
       Temperature = T,
       Treatment = str_replace(Treatment,regex("^evo"),"Evo"))

###
#Read growth from respiration experiments (Flow cytometry)
###
#Read in OTU growth
OTU_r_flow <- read_csv("../data/Francisca/FlowCytometry_All.csv",col_types = cols()) %>%
    filter(Level == "Species", Id != 9) %>%
    mutate(r = (log(AbuT1) - log(AbuT0)) / AT,
           Treatment = ifelse(Treatment == "Denovo","Nonevolved",Treatment),
          OTU = Id) %>%
    select(OTU,Temperature,Treatment,Media,Replicate,r) %>%
    filter(r > -1)

#Read in Pairwise growth
pairs_r_flow <- read_csv("../data/Francisca/FlowCytometry_All.csv", col_types = cols()) %>%
    filter(Level == "Pairs") %>%
    mutate(r = (log(AbuT1)-log(AbuT0))/AT,
           Treatment = recode(Treatment, "Denovo" = "Nonevolved"),
           Temperature = Temperature) %>%
    separate(Id,c("OTU_1","OTU_2"),sep = "-") %>%
    mutate(OTU_1 = ifelse(OTU_1 == "Feb",2,OTU_1)) 

####
#From respiration experiments (OD)
####
#Read in OTU growth
OTU_r_OD_1 <- read_csv("../data/Francisca/respirationandODData/DatosOD_5OTUs_evNev.csv", col_types = cols()) %>%
    mutate(r = (log(OD) - log(OD_T0))/(Tf / 60)) %>%
#            r = ifelse(r < 0 , 0 , r)) %>%
    select(OTU,Temperature = T, Treatment, Replicate, r)

#Read in Pairwise growth
pairs_r_OD <- read_csv("../data/Francisca/respirationandODData/DatosOD_PairsevNev_u.csv", col_types = cols()) %>%
    mutate(r = (log(OD) - log(OD_T0))/(Tf / 60),
#            r = ifelse(r < 0 , 0 , r),
           OTU = str_replace(OTU,"Feb","2"),
           Treatment = ifelse(Treatment == "evolved", "Evolved",Treatment)) %>%
    separate(OTU,into = c("OTU_1","OTU_2"),sep = "-") %>%
    select(OTU_1,OTU_2,Temperature = T, Treatment, Replicate, r,AbuT0 = OD_T0)
    
####
#From new full growth curve experiments
####
OTU_r_OD_2 <- read_csv("../data/Francisca/LogisticIndivcurves_out_M9_updated.csv", 
                     col_types = cols()) %>%
#     filter(!is.na(Q_r)) %>%
    separate(pa,c("OTU","Media","Treat"),"-") %>%
    separate(Treat,c("Treat","Temp","Rep"),":") %>%
    mutate(Temperature = as.numeric(Temp),
           Temperature = ifelse(Temperature == 27, 27.5, Temperature),
           Treatment = recode(Treat, "NE" = "Nonevolved", "E" = "Evolved"),
           OTU = str_remove(OTU,"OTU")) %>%
#            r = ifelse(r < 0 , 0 , r)) %>%
    select(OTU,Media,Temperature,Treatment,r,K)


###
# K (carrying capacity) @ 20 degrees
###
#getting K estimates
K_data <- read_csv("../data/Francisca/GrowthRate20C/Pairwise_LogisticIndivcurves_Barout.csv",
         col_types = cols()) %>%
    filter(str_detect(Treatment,"\\.",negate = TRUE)) %>%
    separate(Treatment, into = c("OTU","Treatment"), sep = "-") %>%
    mutate(K = 10^(LOG10Nmax)) %>%
    filter(K < 0.75) %>%
    group_by(OTU,Treatment) %>%
    summarise(uK = mean(K), .groups = "drop") %>%
    ungroup() %>%
    mutate(Treatment = recode(Treatment,E = "Evolved", NE = "Nonevolved"))


#------------------------------------------------------------

###
# plot data distributions
###

cPalette <- c("#F0E442", "#0072B2", "#D55E00")[c(3,2)]

###
#Flow cytometry data
###
OTU_flow = OTU_r_flow %>% select(Treatment,r)
pair_flow = pairs_r_flow %>% select(Treatment,r)
Data_source = c(rep("OTU_flow",nrow(OTU_flow)), rep("pair_flow",nrow(pair_flow)))

p6 = bind_rows(OTU_flow,pair_flow) %>%
    mutate(Data_source) %>%
    ggplot(aes(x=r,fill = Treatment))+
         geom_density(position = "identity",alpha = 0.7)+
         geom_vline(xintercept = 0.0)+
         scale_fill_manual(values = cPalette)+   
         theme_cowplot() +
         theme(legend.position = "none",plot.margin = margin(25,0,0,0))+
         facet_wrap(~Data_source,ncol = 1)



###
#OD data
###
OTU_OD_1 = OTU_r_OD_1 %>% select(Treatment,r)
OTU_OD_2 = OTU_r_OD_2 %>% select(Treatment,r)
pair_OD = pair_flow %>% select(Treatment,r)
Data_source = c(rep("OTU_OD_1",nrow(OTU_OD_1)), rep("OTU_OD_2",nrow(OTU_OD_2)), rep("pair_OD",nrow(pair_OD)))

p7 <- bind_rows(OTU_OD_1,OTU_OD_2,pair_OD) %>%
    mutate(Data_source) %>%
    filter(r > -2.5 , r < 2.5) %>%
    ggplot(aes(x= r ,fill = Treatment))+
        geom_density(position = "identity",alpha = 0.7)+
        geom_vline(xintercept = 0.0)+
        scale_fill_manual(values = cPalette)+
        theme_cowplot()+
        theme(legend.title = element_blank(),,plot.margin = margin(25,0,0,0))+
        facet_wrap(~Data_source,ncol = 1)

#----------------------------------------------------------
###
#plot data over temperature
###




#----------------------------------------------------------

###
# Do the bootstraping to get interaction estimates
###

#Flow cytometry interaction estimates

#get unique combinations
Temp_Treat_df = pairs_r_flow %>%
    group_by(Temperature,Treatment) %>%
    summarise(.groups = 'drop') 

grw_list_flow <- list()
res_list_flow <- list()
int_list_flow <- list()

for(i in 1:nrow(Temp_Treat_df)){
    estimates = bootstrap(Temp_Treat_df$Treatment[i],Temp_Treat_df$Temperature[i],100, "flow","flow",K_data)
    grw_list_flow[[i]] = estimates[[1]]
    res_list_flow[[i]] = estimates[[2]]
    int_list_flow[[i]] = estimates[[3]]
}

#OD interaction estimates
#get unique combinations
Temp_Treat_df_2 <- OTU_r_OD_2 %>% 
    group_by(Temperature,Treatment) %>%
    summarise(.groups = 'drop') %>%
    ungroup() %>%
    semi_join(Temp_Treat_df)

grw_list_OD <- list()
res_list_OD <- list()
int_list_OD <- list()

for(i in 1:nrow(Temp_Treat_df_2)){
    estimates = bootstrap(Temp_Treat_df_2$Treatment[i],Temp_Treat_df_2$Temperature[i],1000, "OD_1","OD_1",K_data)
    grw_list_OD[[i]] = estimates[[1]]
    res_list_OD[[i]] = estimates[[2]]
    int_list_OD[[i]] = estimates[[3]]
}
