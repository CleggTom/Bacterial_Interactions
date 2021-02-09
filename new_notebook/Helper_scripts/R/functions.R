#### FUNCTIONS.R
# This script contains helper functions for the supplimentary notebook

###
#Functions to do bootstraping
###
#function to get bootstrap samples for OTU respiration rates
get_OTU_R_samples = function(OTU_vec, Treat, Temp, n){
    #get individual samples
    OTU_samples <- list()
    for(i in 1:length(OTU_vec)){
        OTU_samples[[i]] = OTU_R %>% 
            filter(OTU == OTU_vec[i], Treatment == Treat, Temperature == Temp) 
        
        if(nrow(OTU_samples[[i]]) > 0){
          OTU_samples[[i]] = OTU_samples[[i]] %>%
                pull(Resp) %>%
                sample(.,n,TRUE)
        } else {
            OTU_samples[[i]] = rep(0,n) 
        }
    }
    names(OTU_samples) <- OTU_vec 
    return(OTU_samples)
}

#function to get a bootstraped sample for single OTU growth rates
get_OTU_r_samples = function(OTU_vec, Treat, Temp, n, data = ""){
    #choose dataset 
    if(data == "flow"){
     x = OTU_r_flow  
    } else if(data == "OD_1") {
     x = OTU_r_OD_1   
    } else if(data == "OD_2"){
     x = OTU_r_OD_2  
    } else {stop("dataset not specified correctly")}
        
    #get individual samples
    OTU_samples <- list()
    for(i in 1:length(OTU_vec)){
        OTU_samples[[i]] = x %>% 
            filter(OTU == OTU_vec[i], Treatment == Treat, Temperature == Temp) %>%
            pull(r) %>%
            sample(.,n,TRUE)
    }
    names(OTU_samples) <- OTU_vec 
    return(OTU_samples)
}

get_OTU_r_T = function(Treat, data, T){
    if(data == "flow"){
     x = OTU_r_flow  
    } else if(data == "OD_1") {
     x = OTU_r_OD_1   
    } else if(data == "OD_2"){
     x = OTU_r_OD_2  
    } else {stop("dataset not specified correctly")}
        
    x %>%
        filter(Temperature == T,Treatment == Treat) %>%
        group_by(OTU) %>%
        summarise(mean_r = mean(r),.groups = 'drop') %>%
        spread(OTU,mean_r)  
}

#function to get bootstraped paired growth and biomasss
get_Paris_r_samples = function(OTU1,OTU2,Treat, Temp, n,data = ""){
    #choose dataset 
    if(data == "flow"){
     x = pairs_r_flow 
    } else if(data == "OD_1") {
     x = pairs_r_OD   
    } else {stop("dataset not specified correctly")}
        
    #get pair samples
    Pairs_r_samples <- list()
    Pairs_N_samples <- list()

    for(i in 1:length(OTU1)){
        y = x %>%
            filter(OTU_1 == OTU1[i], OTU_2 == OTU2[i], Treatment == Treat, Temperature == Temp)
        
        indx = sample(1:nrow(y),n,TRUE)
        Pairs_r_samples[[i]] = y$r[indx]
        Pairs_N_samples[[i]] = y$AbuT0[indx]
    }
    names(Pairs_r_samples) <- paste(OTU1,OTU2,sep = "-")
    names(Pairs_N_samples) <- paste(OTU1,OTU2,sep = "-")
    
    return(list(Pairs_r_samples,Pairs_N_samples))
}

#function to generate bootstrapped communities at a given temp/treatment
bootstrap <- function(Treat,Temp,n, OTU_data, pair_data, K_data){
    #bootstrap OTU R
    OTU_vec = unique(OTU_R$OTU)
    OTU_R_sample = get_OTU_R_samples(OTU_vec,Treat,Temp,n)
    
    #bootstrap OTU r
    OTU_vec = unique(OTU_r_flow$OTU)
    OTU_r_sample = get_OTU_r_samples(OTU_vec,Treat,Temp,n,OTU_data)
    
    #get individual growth at 20C
    OTU_r_T = get_OTU_r_T(Treat,OTU_data,15)

    #boostrap pairs r
    iter_df = pairs_r_flow %>%
        group_by(OTU_1,OTU_2) %>%
        summarise(.groups = 'drop') 
    Pairs_sample = get_Paris_r_samples(iter_df$OTU_1,iter_df$OTU_2,Treat,Temp,n,pair_data)

    #calculate a values
    a_list <- list()
    #loop through all pairs
    for(i in 1:nrow(iter_df)){
        r1 = OTU_r_sample[[which(iter_df$OTU_1[i] == OTU_vec)]]
        r2 = OTU_r_sample[[which(iter_df$OTU_2[i] == OTU_vec)]]
        
        rp = Pairs_sample[[1]][[i]]
        cp = Pairs_sample[[2]][[i]]

        K1 = K_data %>% filter(OTU == iter_df$OTU_1[i], Treatment == Treat) %>%
                select(uK) %>% unlist %>% unname
        K2 = K_data %>% filter(OTU == iter_df$OTU_2[i], Treatment == Treat) %>%
                select(uK) %>% unlist %>% unname
        
        a_list[[i]] = ( ((2*rp) - r1 - r2) / cp) - 0.0

    }
    names(a_list) <- names(Pairs_sample[[1]])

    #combine to dataframes
    #growth rates 
    growth_df = bind_cols(OTU_r_sample) %>%
        mutate(Rep = 1:n) %>%
        gather("OTU","r",-Rep) %>%
        mutate(Temperature = Temp, Treatment = Treat)

    #resp
    resp_df = bind_cols(OTU_R_sample) %>%
        mutate(Rep = 1:n) %>%
        gather("OTU","R",-Rep) %>%
        mutate(Temperature = Temp, Treatment = Treat)    
    
    #get pairs    
    int_df = bind_cols(a_list) %>%
        mutate(Rep = 1:n) %>%
        gather("pair","a", -Rep) %>%
        separate(pair,c("OTU_1","OTU_2"),"-")%>%
        mutate(Temperature = Temp, Treatment = Treat)
    
    return(list(growth_df,resp_df,int_df))
}