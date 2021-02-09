#### SIMULATION ANALYSIS.R
# This script contains the code to do analysis and plotting of data from the notebook

#--------------------------------------------------

###
#Plotting simulated data
###

#function to get intervals of continuous columns
sampled_mode <- function(vec,p){
    x <- sort(unique(vec))
    indx <- ceiling(p * length(x))
    return(x[indx])
}

Cfiles <- list.files("../Data/Simulations/Theory/C",full.names = T)
Rfiles <- list.files("../Data/Simulations/Theory/R",full.names = T)
afiles <- list.files("../Data/Simulations/Theory/a",full.names = T)

Nint = 50
NE = 10
Ntemp = 15
Ntime = 50

T <- seq(1,-1,length.out = Ntemp)
t <- seq(0,10,length.out = Ntime)
a <- seq(-3,1,length.out = Nint)
Evec <- seq(0,2,length.out = NE)

df_list <- list()

for(i in 1:length(Rfiles)){
    x <- suppressMessages(read_csv(Rfiles[i]))
    y <- suppressMessages(read_csv(afiles[i]))
    z <- suppressMessages(read_csv(Cfiles[i]))
    
    colnames(x) <- T
    colnames(y) <- T
    colnames(z) <- T

    #get Respiration df
    x <- x %>% mutate(E = Evec[i],
                      Time = rep(t, each = Nint),
                      tN = rep(1:Ntime, each = Nint),
                      aN = rep(1:Nint, times = Ntime)) %>%
           gather("Temp","R",-E,-Time,-aN,-tN) %>%
            mutate(Temp = as.numeric(Temp)) 
    
    #get interaction df
    y <- y %>% mutate(E = Evec[i],  Time = rep(t, each = Nint)) %>%
           gather("Temp","a",-E,-Time) %>%
            mutate(Temp = as.numeric(Temp)) 
    
    #get biomass df
    z <- z %>% mutate(E = Evec[i],  Time = rep(t, each = Nint)) %>%
           gather("Temp","C",-E,-Time) %>%
            mutate(Temp = as.numeric(Temp)) 

    x$a <- y$a
    x$C <- z$C
    df_list[[i]] <- x
}

####
#interactions
###

p1 <- bind_rows(df_list) %>%
    filter(E == 0.0, Temp == 0) %>%
    ggplot(aes(x=Time,y=C,colour=a, group = aN))+
            geom_line() + 
            geom_line(data =  bind_rows(df_list) %>% 
                                  filter(E == 0.0, Temp == 0, aN == 25),
                      colour = "red" , size = 2)+
            theme_cowplot()

p2 <- bind_rows(df_list) %>%
    filter(E == 0.0, Temp == 0) %>%
    ggplot(aes(x=Time,y=R,colour=a, group = aN))+
            geom_line() + 
            geom_line(data =  bind_rows(df_list) %>% 
                                  filter(E == 0.0, Temp == 0, aN == 25),
                      colour = "red" , size = 2)+
            theme_cowplot()
###
#temperature
###
p3 <- bind_rows(df_list) %>%
    filter(E == 0 , aN == 25) %>%
    ggplot(aes(x=Time, y = C, group = Temp, colour = Temp))+
        geom_line()+
        geom_vline(xintercept = 1)+
        scale_color_gradient(low="red", high="blue",guide = F)+
        theme_cowplot()

p4 <- bind_rows(df_list) %>%
    filter(E == 0 ,aN == 25 , tN == 10 ) %>%
    ggplot(aes(x=-Temp,y=C))+
        geom_point()+
        theme_cowplot()

###
# Across Interactions strengths
###
p5 <- bind_rows(df_list) %>%
    filter(E == 0.0,tN == 10) %>%
    ggplot(aes(x=Temp, y = (R), colour = a))+
        geom_point() +
        theme_cowplot()
