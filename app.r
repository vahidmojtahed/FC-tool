## this is a single Shiny file version - although we keep model.r as a separate file for convenience ----

## install if needed and load required libraries ----
if (!"deSolve" %in% installed.packages()) install.packages("deSolve")
library(deSolve)
if (!"shiny" %in% installed.packages()) install.packages("shiny")
library(shiny)
if (!"shinythemes" %in% installed.packages()) install.packages("shinythemes")
library(shinythemes)
if (!"triangle" %in% installed.packages()) install.packages("triangle")
library(triangle)

## load the model ----
source("model.r")

## initiate the random generators ----
globalSeed<-NULL

## helper functions for plotting ----
density.new<-function (data, nclass = 21, probability = F,...) 
{
  xx <- pretty(data, n = nclass)
  xx.h <- hist(data, breaks = xx, plot = F)
  xx.h.x <- xx.h$breaks
  xx.h.x <- xx.h.x[1:(length(xx.h.x) - 1)] + 0.5 * diff(xx.h.x)
  if (probability) {
    xx.h.y <- xx.h$density
  }
  else {
    xx.h.y <- xx.h$counts
  }
  return(list(x = xx.h.x, y = xx.h.y))
}

## for setting and saving random seed - ensuring that the simulations are random but reproducible ----
# the implementation is not as generic as I would like to, as it precludes stochastic models
# this is done whenever variability() is called and a new set of params.run is generated
# the seed is saved in parameters - currently not used
seed.init<-function(globalSeed){
  if (is.null(globalSeed)){ 
    as.numeric(Sys.time())-> t; set.seed((t - floor(t)) * 1e8 -> seed) 
    } else seed<-globalSeed
  return(seed)
}

## Server - this is where all calculations take place ----
server<-function(input, output) {
  
  ## base parameters ----
  baseParams<-reactive({
    # fixed through scenario
    params<-selectParams()
    N<-params[1]
    beta.base<-params[2] # to save the current value, we need to read in beta not beta base - not 4
    y0.base<-params[3] # ditto for y0 - not 8
#    value<-params[4]
#    contr<-params[5] # effectiveness of control
    value.timber<-params[4]
    value.recreation<-params[5]
    value.landscape<-params[6]
    value.biodiversity<-params[7]
    value.carbon<-params[8]
    contr<-params[9] # effectiveness of control
    control.costs<-params[10]
    Nrep<-params[13]
    d.rate<-params[14]
    globalSeed<-seed.init(globalSeed)
    return(list(
      N=N,
      control.costs=control.costs,
      Nrep=Nrep,
      d.rate=d.rate,
      seed=globalSeed,
      beta.base=beta.base,
      y0.base=y0.base,
      contr=contr,
      value.timber=value.timber,
      value.recreation=value.recreation,
      value.landscape=value.landscape,
      value.biodiversity=value.biodiversity,
      value.carbon=value.carbon
    ))
  })

check.in<-function(xin,xbase){
  return(ifelse(is.null(xin),xbase,xin))
#  return(xbase)
}
  
## load parameters and update them in real time ----     
  get.params<-function(){
    params.base<-baseParams()
    Nrep<-check.in(input$Nrep,params.base$Nrep)
    N<-check.in(input$Nsize,params.base$N)
    d.rate<-check.in(input$d.rate,params.base$d.rate)
    if (is.null(input$control.varmodel)) control.varmodel<-1 else control.varmodel<-input$control.varmodel
    if (is.null(input$control.prop)) control.prop<-0 else control.prop<-input$control.prop
    return(list(
      N=N,        
      beta=params.base$beta.base*mean(input$beta.var)/100, # for triangular assume it is symmetric
      y0=params.base$y0.base*mean(input$y0.var)/100,
      value.timber=check.in(input$Timber,params.base$value.timber),
      value.recreation=check.in(input$Recreation,params.base$value.recreation),
      value.landscape=check.in(input$Landscape,params.base$value.landscape),
      value.biodiversity=check.in(input$Biodiversity,params.base$value.biodiversity),
      value.carbon=check.in(input$Carbon,params.base$value.carbon),
      contr=params.base$contr,
      control.cost=check.in(input$Cost,params.base$control.cost),
      control.varmodel=control.varmodel,
      control.prop=control.prop,
      Nrep=Nrep,
      d.rate=d.rate,
      effort=input$effort,
      beta.base=params.base$beta.base,
      beta.var1=params.base$beta.base*input$beta.var[1]/100,
      beta.var2=params.base$beta.base*input$beta.var[2]/100,
      y0.base=params.base$y0.base,
      y0.var1=params.base$y0.base*input$y0.var[1]/100,
      y0.var2=params.base$y0.base*input$y0.var[2]/100,
      seed=params.base$seed
    ))
  }
  
  ## main simulation engine  ----
  DEoutput<-function(params){
    # solve equations
    # this is very badly coded
#    params<-get.params()
    Nrep1<-1:params$Nrep             # replicates with control
    Nrep2<-params$Nrep+1:params$Nrep        # replicates without control
    set.seed(params$seed)
    S<-matrix(NA,nrow=length(times),ncol=2*params$Nrep)
    I<-matrix(NA,nrow=length(times),ncol=2*params$Nrep)
    for (i in Nrep1){
      params.run<-variability(params,contr.switch=1)
      yinit<-c(S=params$N-params.run$y0,I=params.run$y0)
      result<-ode(yinit,times,forest,params.run,method="bdf")
      S[,i]<-result[,2]
      I[,i]<-result[,3]
    }
    for (i in Nrep2){
      params.run<-variability(params,contr.switch=0)
      yinit<-c(S=params$N-params.run$y0,I=params.run$y0)
      result<-ode(yinit,times,forest,params.run,method="bdf")
      S[,i]<-result[,2]
      I[,i]<-result[,3]
    }
    return(list(S=S,I=I))
  }
  
  ## runtime; so that we can choose the type of the model - at the moment this simply calls DEoutput ----
  get.output<-function(params){
    return(DEoutput(params))
  }

  ## inputs ----
  output$selectNsize <- renderUI({
    params.base<-get.params()
    numericInput("Nsize","Area (ha):",params.base$N,min=min(params.base$N,50),max=max(params.base$N,1000),step=50)
  })
  output$selectNrep <- renderUI({
    params.base<-baseParams()
    numericInput("Nrep","Number of replicates:",params.base$Nrep,min=50,max=1000,step=50)
  })
  output$selectCost <- renderUI({
    params.base<-baseParams()
    numericInput("Cost","Cost per man-week:",params.base$control.costs,min=0,max=100,step=0.1)
  })
  output$selectd.rate <- renderUI({
    params.base<-baseParams()
    numericInput("d.rate","Discount rate (per year):",params.base$d.rate,min=0,max=0.5,step=0.01)
  })
  output$selectTimber <- renderUI({
    params.base<-baseParams()
    numericInput("Timber","Timber (GBP per ha):",params.base$value.timber,min=50,max=10000,step=50)
  })
  output$selectRecreation <- renderUI({
    params.base<-baseParams()
    numericInput("Recreation","Recreation (GBP per ha per yr):",params.base$value.recreation,min=0,max=100,step=5)
  })
  output$selectLandscape <- renderUI({
    params.base<-baseParams()
    numericInput("Landscape","Landscape (GBP per ha per yr):",params.base$value.landscape,min=0,max=100,step=5)
  })
  output$selectBiodiversity <- renderUI({
    params.base<-baseParams()
    numericInput("Biodiversity","Biodiversity (GBP per ha per yr):",params.base$value.biodiversity,min=0,max=100,step=5)
  })
  output$selectCarbon <- renderUI({
    params.base<-baseParams()
    numericInput("Carbon","Carbon (GBP per ha per yr):",params.base$value.carbon,min=0,max=100,step=5)
  })

  
  ## outputs ----
  # traces ----
  output$tracePlot <- renderPlot({
    params<-get.params()
    Nrep0<-1:min(10,params$Nrep)
    out<-get.output(params)
    par(mfrow = c(2, 1), mar = c(3, 1, 1, 1), oma = c(0, 3,0, 0))
    par(tcl = -0.2, pty = "m", cex = 1, las = 1, cex.axis = 1, 
        cex.lab = 1, mgp = c(3, 0.5, 0))
    plot(times,out$S[,1],type="n",ylim=c(0,params$N),xlab="",ylab="Healthy/Infected trees (%)")
    matlines(times,out$I[,Nrep0],type="l",col=2,lty=1)
    matlines(times,out$S[,Nrep0],type="l",col=3,lty=1)
    abline(h=0,lty=2,col=4)
    if (params$effort>0) abline(v=time.stop,lty=2,col=4)
    z<-damage(out$S[,Nrep0],out$I[,Nrep0],params,contr.switch=1,control.prop=params$control.prop)
    matplot(times,z,ylim=range(0,z),type="l",xlab="",ylab="damage (weekly)",lty=1,col=1)
    abline(h=0,lty=2,col=4)
    if (params$effort>0) abline(v=time.stop,lty=2,col=4)
    mtext(side=1,line=2,"Time (years)",cex=1)
  })
  
  # distributions of parameters ----
  output$parsPlot <- renderPlot({
    params<-get.params()
    params.control<-variability(params,contr.switch=1,runs=1000)
    #    params.no<-variability(params,contr.switch=0,runs=500)
    par(mfrow = c(3, 1), mar = c(3, 1, 1, 1), oma = c(0, 3,0, 0))
    par(tcl = -0.2, pty = "m", cex = 1, las = 1, cex.axis = 1, 
        cex.lab = 1, mgp = c(2, 0.5, 0))
    z<-density.new(params.control$beta,adjust=1)
    plot(z$x,z$y/max(z$y),xlim=range(0,z$x),col=1,lwd=2,type="s",ylab="Density",xlab="rate of spread")
    z<-density.new(params.control$y0,adjust=0.5)
    plot(z$x,z$y/max(z$y),xlim=range(0,params.control$N,z$x),col=1,lwd=2,type="s",ylab="Density",xlab="initial load")
    z<-density.new(params.control$effort*params.control$contr,adjust=1)
    plot(z$x,z$y/max(z$y),xlim=range(0,max(z$x),0.5),col=1,lwd=2,type="s",ylab="Density",xlab="control rate")
  })
  
  # distribution of damage ----
  output$distPlot <- renderPlot({
    params<-get.params()
    Nrep1<-1:params$Nrep             # replicates with control
    Nrep2<-params$Nrep+1:params$Nrep        # replicates without control
    out<-get.output(params)
    par(mfrow = c(1, 1), mar = c(3, 1, 1, 1), oma = c(0, 3,0, 0))
    par(tcl = -0.2, pty = "m", cex = 1, las = 1, cex.axis = 1, 
        cex.lab = 1, mgp = c(3, 0.5, 0))
    z<-density(damage(out$S[,Nrep1],out$I[,Nrep1],params,1,control.prop=params$control.prop)[time.end,],adjust=1) # was 0.5
    z0<-density(damage(out$S[,Nrep2],out$I[,Nrep2],params,0,control.prop=params$control.prop)[time.end,],adjust=1)
    plot(z$x,z$y/max(z$y),xlim=range(z$x,z0$x),type="s",main="",col=2,lwd=3)
    lines(z0$x,z0$y/max(z0$y),type="s",col=4,lwd=3)
  })
  
  # a copy for report ----
  output$distPlot.report <- renderPlot({
    params<-get.params()
    Nrep1<-1:params$Nrep             # replicates with control
    Nrep2<-params$Nrep+1:params$Nrep        # replicates without control
    out<-get.output(params)
    par(mfrow = c(1, 1), mar = c(3, 1, 1, 2), oma = c(0, 1,0, 1))
    par(tcl = -0.2, pty = "m", cex = 1, las = 1, cex.axis = 1, 
        cex.lab = 1, mgp = c(3, 0.5, 0))
    z<-density(damage(out$S[,Nrep1],out$I[,Nrep1],params,1,control.prop=params$control.prop)[time.end,],adjust=1) # was 0.5
    z0<-density(damage(out$S[,Nrep2],out$I[,Nrep2],params,0,control.prop=params$control.prop)[time.end,],adjust=1)
    plot(z$x,z$y/max(z$y),xlim=range(z$x,z0$x),type="s",main="",col=2,lwd=3)
    lines(z0$x,z0$y/max(z0$y),type="s",col=4,lwd=3)
  })
  
  # table ----
  output$table <- renderTable({
    params<-get.params()
    Nrep1<-1:params$Nrep             # replicates with control
    Nrep2<-params$Nrep+1:params$Nrep        # replicates without control
    out<-get.output(params)
    x<-mean(out$S[time.end,Nrep1],na.rm=T)
    y<-mean(out$I[time.end,Nrep1],na.rm=T)
    y0<-mean(out$I[time.end,Nrep2],na.rm=T)
    z<-mean(damage(out$S[,Nrep1],out$I[,Nrep1],params,1,control.prop=params$control.prop)[time.end,])
    z0<-mean(damage(out$S[,Nrep2],out$I[,Nrep2],params,0,control.prop=params$control.prop)[time.end,])
    x<-rbind(      
      c("Healthy stock (%):",round(100*x/params$N,0)),
      c("Infected stock (%):",round(100*y/params$N,0)),
      c("damage (control):",round(z,0)),
      c("damage (do nothing):",round(z0,0))
    )
    colnames(x)<-c("Indicator","Value")
    return(x)
  },bordered=TRUE)
  
  # table for report ----
  output$table.report <- renderTable({
    params<-get.params()
    Nrep1<-1:params$Nrep             # replicates with control
    Nrep2<-params$Nrep+1:params$Nrep        # replicates without control
    out<-get.output(params)
    x<-mean(out$S[time.end,Nrep1],na.rm=T)
    y<-mean(out$I[time.end,Nrep1],na.rm=T)
    y0<-mean(out$I[time.end,Nrep2],na.rm=T)
    z<-mean(damage(out$S[,Nrep1],out$I[,Nrep1],params,1,control.prop=params$control.prop)[time.end,])
    z0<-mean(damage(out$S[,Nrep2],out$I[,Nrep2],params,0,control.prop=params$control.prop)[time.end,])
    x<-rbind(      
      c("Healthy stock (%):",round(100*x,0)),
      c("Infected stock (%):",round(100*y,0)),
      c("damage (control):",round(z,0)),
      c("damage (do nothing):",round(z0,0))
    )
    colnames(x)<-c("Indicator","Value")
    return(x)
  },bordered=TRUE)
  
  # table showing (some) parameters ----
  output$parametersTable <- renderTable({
    params<-get.params()
    x<-rbind(
      c("Total area:",round(params$N,2)),
      c("Rate of spread:",round(params$beta.base,2)),
      c("Initial infected population:",round(params$y0,2)),
      c("Cost of control:",round(params$control.cost,5))
    )
    colnames(x)<-c("Parameter","Value")
    return(x)
  },bordered=TRUE)
  
  # same for report ----
  output$parametersTable.report <- renderTable({
    params<-get.params()
    x<-rbind(
      c("Total area:",round(params$N,2)),
      c("Rate of spread:",round(params$beta.base,2)),
      c("Initial infected population:",round(params$y0,2)),
      c("Cost of control:",round(params$control.cost,5))
    )
    colnames(x)<-c("Parameter","Value")
    return(x)
  },bordered=TRUE)
  
  # show scenarios for upload tab ----
  output$currentscenarioTable.scenarios <- renderTable({
    params<-get.params()
      x<-data.frame("ID"=1,
                    "Pest"="Generic",
                    "Host"="Generic",
                    "UK"="Present",
                    "Regulation"="FC",
                    "Likelihood"=1,
                    "Impact"=3,
                    "Value"=1,
                    "Risk"=40,
                    "N"=params$N,
                    "beta"=params$beta,
                    "y0"=params$y0,
                    "value.timber"=params$value.timber,
                    "value.recreation"=params$value.recreation,
                    "value.landscape"=params$value.landscape,
                    "value.biodiversity"=params$value.biodiversity,
                    "value.carbon"=params$value.carbon,
                    "contr"=params$contr,
                    "control.cost"=params$control.cost,
                    "control.varmodel"=params$control.varmodel,
                    "control.prop"=params$control.prop,
                    "Nrep"=params$Nrep,
                    "d.rate"=params$d.rate,
                    "effort"=params$effort,
                    "beta.base"=params$beta.base,
                    "beta.var1"=params$beta.var1,
                    "beta.var2"=params$beta.var2,
                    "y0.base"=params$y0.base,
                    "y0.var1"=params$y0.var1,
                    "y0.var2"=params$y0.var2,
                    "seed"=params$seed
      )[10:31]
      colnames(x)<-c("ID","Pest","Host","UK","Regulation","Likelihood","Impact","Value","Risk",
                     "N","beta","y0","value.timber","value.recreation","value.landscape","value.biodiversity","value.carbon","contr","control.cost","control.varmodel","control.prop","Nrep","d.rate","effort",
                     "beta.base","beta.var1","beta.var2","y0.base","y0.var1","y0.var2","seed")[10:31]
      return(x)
  },bordered=TRUE)
  output$currentscenarioTable.upload <- renderTable({
    params<-get.params()
    x<-data.frame("ID"=1,
                  "Pest"="Generic",
                  "Host"="Generic",
                  "UK"="Present",
                  "Regulation"="FC",
                  "Likelihood"=1,
                  "Impact"=3,
                  "Value"=1,
                  "Risk"=40,
                  "N"=params$N,
                  "beta"=params$beta,
                  "y0"=params$y0,
                  "value.timber"=params$value.timber,
                  "value.recreation"=params$value.recreation,
                  "value.landscape"=params$value.landscape,
                  "value.biodiversity"=params$value.biodiversity,
                  "value.carbon"=params$value.carbon,
                  "contr"=params$contr,
                  "control.cost"=params$control.cost,
                  "control.varmodel"=params$control.varmodel,
                  "control.prop"=params$control.prop,
                  "Nrep"=params$Nrep,
                  "d.rate"=params$d.rate,
                  "effort"=params$effort,
                  "beta.base"=params$beta.base,
                  "beta.var1"=params$beta.var1,
                  "beta.var2"=params$beta.var2,
                  "y0.base"=params$y0.base,
                  "y0.var1"=params$y0.var1,
                  "y0.var2"=params$y0.var2,
                  "seed"=params$seed
    )[10:31]
    colnames(x)<-c("ID","Pest","Host","UK","Regulation","Likelihood","Impact","Value","Risk",
                   "N","beta","y0","value.timber","value.recreation","value.landscape","value.biodiversity","value.carbon","contr","control.cost","control.varmodel","control.prop","Nrep","d.rate","effort",
                   "beta.base","beta.var1","beta.var2","y0.base","y0.var1","y0.var2","seed")[10:31]
    return(x)
  },bordered=TRUE)
  
    output$scenarioTable.upload <- renderTable({
    params<-get.params()
    pars.in<-parsInput()[2:31]
    if (is.null(pars.in)) 
      return() 
    else {
      x<-pars.in
      colnames(x)<-c("ID","Pest","Host","UK","Regulation","Likelihood","Impact","Value","Risk",
                     "N","beta","y0","value.timber","value.recreation","value.landscape","value.biodiversity","value.carbon","contr","control.cost","control.varmodel","control.prop","Nrep","d.rate","effort",
                     "beta.base","beta.var1","beta.var2","y0.base","y0.var1","y0.var2","seed")[2:31]
      return(x)
    }
  },bordered=TRUE)
  
  # scenarios for scenarios tab ----
  output$scenarioTable.scenarios <- renderTable({
    params<-get.params()
    pars.in<-parsInput()[2:31]
    if (is.null(pars.in)) 
      return() 
    else {
      x<-pars.in
      colnames(x)<-c("ID","Pest","Host","UK","Regulation","Likelihood","Impact","Value","Risk",
                     "N","beta","y0","value.timber","value.recreation","value.landscape","value.biodiversity","value.carbon","contr","control.cost","control.varmodel","control.prop","Nrep","d.rate","effort",
                     "beta.base","beta.var1","beta.var2","y0.base","y0.var1","y0.var2","seed")[2:31]
      return(x)
    }
  },bordered=TRUE)
  
  ## helpers to output data ----
  # dump simulation results ----
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$dataset, '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file)
    }
  )
  datasetInput <- reactive({
    params<-get.params()
    Nrep1<-1:params$Nrep             # replicates with control
    Nrep2<-params$Nrep+1:params$Nrep        # replicates without control
    out<-get.output(params)
    switch(input$dataset,
           "Healthy" = out$S[,1+Nrep1],
           "Infected" = out$I[,1+Nrep1])
  })
  
  # dump parameters/scenarios ----
  output$downloadPars <- downloadHandler(
    filename = function() { paste("State", '.csv', sep='') },
    content = function(file) {
      write.csv(paramsInput(), file)
    }
  )
  paramsInput <- reactive({
    params<-get.params()
    return(c(selectParamsNames(),params))
  })
  
  ## management of scenarios ----
  # upload scenarios ----
  parsInput <- reactive({
    inFile <- input$uploadPars
    if (is.null(inFile)){
      uploadPars<-NULL
    } else {
      uploadPars<-read.csv(inFile$datapath,header=F,skip=1,                           ,
                             col.names=c(
                               "ID","Pest","Host","UK","Regulation","Likelihood","Impact","Value","Risk",
                               "N","beta","y0","value.timber","value.recreation","value.landscape","value.biodiversity","value.carbon","contr","control.cost","control.varmodel","control.prop","Nrep","d.rate","effort",
                               "beta.base","beta.var1","beta.var2","y0.base","y0.var1","y0.var2","seed"
                            )
                  )
    }
    return(rbind.data.frame(defaultScenario(), uploadPars))
  })
  

  # select scenarios ----
  output$selectScenario <- renderUI({
    scores.in<-parsInput()
    scores.choices<-setNames(1:nrow(scores.in),scores.in[,2])
    selectInput("scenarios", 
                label = h5("Select scenario:"), 
                choices = scores.choices,
                selected = 1)
  })
  
  # populate parameters ----
  # this really should be done better, based on column names not on hard wired column numbers
  selectParams <- reactive({
    scores.in<-parsInput()
    if (is.null(input$scenarios)) 
      scores.base<-as.numeric(scores.in[1,10:31])
    else
      scores.base<-as.numeric(scores.in[input$scenarios,10:31])
    return(scores.base)
  })
  
  # populate the rest of the scenarios to match Plant Health Risk Register format ----
  selectParamsNames <- reactive({
    scores.in<-parsInput()
    if (is.null(input$scenarios)) 
      scores.base<-scores.in[1,2:9]
    else
      scores.base<-scores.in[input$scenarios,2:9]
    return(scores.base)
  })
  
  # helper for selecting scenarios ----
  output$test.scenarios <- renderText({
    if (is.null(input$scenarios)) c(1,selectParams())
    else c(paste(input$scenarios),selectParams())
  })
  
}


## User interface part - this is where the output is constructed ----
ui <- navbarPage("FC Tool",id="select.tab",
           theme=shinytheme("superhero"),
           #           shinythemes::themeSelector(),
## welcome tab ----
           tabPanel("Welcome",
                    helpText(h4("Forest disease decision support tool")),
                    helpText("Created by Adam Kleczkowski and Morag Macpherson"),
                    helpText("Project led by Glyn Jones, FERA"),
                    helpText("Collaboration with Julia Touza and Piran White, York, and Stephen Parnell, Salford"),
                    tags$hr(),
                    helpText(h5("Warning: This is a development version. Only limited checking of variables and parameters is performed, hence the tool can break down or produce unreliable results if parameters are entered outside the realistic range.")),
                    tags$hr(),
                    helpText(h5("Version history:")),
                    helpText("Version 0.1: single zone; eradication control; control cost for whole area; value proportional to healthy area."),
                    helpText("Version 0.2: improved GUI and additional inputs/outputs; linked to PHRR"),
                    helpText("Version 0.3: changes to menu; system of equations rather than a single one; discounted economics"),
                    helpText("Version 1.0: moved to a single file; use triangular distributions; changed to years; control only the first 5 years; profit changed to damage; changed from scenarios to parameters"),
                    helpText("Version 1.1: damage function now includes two classes of values; expanded values input; improved front end and saving and loading scenarios")
           ),
## dashboard - where most action takes place ----
           tabPanel("Dashboard",
                    sidebarLayout(
                      sidebarPanel(
                        fluidPage(
                          helpText(h4("Choose the parameters:")),
                          fluidRow(
                            sliderInput("beta.var",
                                        h5("Select min/max for infection rate multiplier (how uncertain are we about future disease spread? 100%=baseline):"),
                                        min = 0,
                                        max = 500,
                                        step=5,
                                        value = c(75,125))
                          ),
                          fluidRow(
                            sliderInput("y0.var",
                                        h5("Select initial state variability (how uncertain are we about how much disease is now? 100%=baseline):"),
                                        min = 1,
                                        max = 500,
                                        step=5,
                                        value = c(75,125))
                          ),
                          fluidRow(
                            sliderInput("effort",
                                        h5("Select control effort (number of man-years to implement control measures?):"),
                                        min = 0,
                                        max = 1000,
                                        step=10,
                                        value = 0.0)
                          ),
                          submitButton("Update")
                        ) # fluid page
                      ), # sidebar panel
                      # Show a plot of the generated distribution
                      mainPanel(
                        fluidPage(
                          fluidRow(
                            #                          helpText(h4("Predictions:")),
                            column(6,
                                   h5("Proportion of healthy and infected trees 
                                      and total damage over a year:"),
                                   h4(tableOutput("table"))
                                   ),
                            column(6,
                                   h5("Distribution of combined damage:\nred=with control, blue=do nothing"),
                                   plotOutput("distPlot",width="275px",height="200px")
                            )
                          ),
                          fluidRow(
                            column(6,
                                   helpText(h4("Simulation trace:")),
                                   #            plotOutput("tracePlot",width="time.end5px",height="350px")
                                   plotOutput("tracePlot",width="350px",height="350px"),
                                   helpText(h5("(green=healthy, red=infected; 10 replicates are shown)"))
                            ),
                            column(6,
                                   helpText(h4("Parameter distributions")),
                                   plotOutput("parsPlot",width="275px",height="350px")
                            )
                          ) # fluid row
                          ) # fluid page
                    ) # main panel
                    )
           ), # Dashboard
## settings menu ----
           navbarMenu("Settings",
                      tabPanel("General parameters",
                               helpText(h4("Set general parameters")),
                               fluidPage(
                                 fluidRow(
                                   column(6,
                                          uiOutput("selectNsize"),
                                          uiOutput("selectCost"),
                                          uiOutput("selectd.rate"),
                                          helpText(h4("Run-time parameters")),
                                          helpText(h5("(set with caution, as large number means slow runs)")),
                                          uiOutput("selectNrep")
                                   ),
                                   column(6,
                                          selectInput("control.prop", 
                                                      label = h5("Select model for control costs:"), 
                                                      choices = c("All area"=0,"Infected area only"=1),
                                                      selected = 0),
                                          selectInput("control.varmodel", 
                                                      label = h5("Select model for control variability:"), 
                                                      choices = c("No variability"=0,"Exponential distribution"=1),
                                                      selected = 0)
                                   )
                                 ) # fluidRow
                               ), # fluidPage
                               tags$hr(),
                               submitButton("Update"),
                               tags$hr(),
                               h4("Current list of key parameters:"),
                               h3(tableOutput("parametersTable"))
                      ),
                      tabPanel("Values at risk",
                               helpText(h4("Set values at risk")),
#                               helpText(h5("Currently implemented via scenarios")),
#                              helpText(h5("At harvest:")),
  fluidRow(
    column(6,
           helpText(h5("Specify the value of timber per ha. This will only be included once per crop (i.e. if trees were harvested at a given time)"))
    ),
    column(6,
          uiOutput("selectTimber")
    )
  ),
fluidRow(
  column(6,
         helpText(h5("Specify the value of recreation:"))
  ),
  column(6,
         uiOutput("selectRecreation")
  )
),
fluidRow(
  column(6,
         helpText(h5("Specify the value of landscape:"))
  ),
  column(6,
         uiOutput("selectLandscape")
  )
),
fluidRow(
  column(6,
         helpText(h5("Specify the value of biodiversity:"))
  ),
  column(6,
         uiOutput("selectBiodiversity")
  )
),
fluidRow(
  column(6,
         helpText(h5("Specify the value of carbon:"))
  ),
  column(6,
         uiOutput("selectCarbon")
  )
),
                              tags$hr(),
                               submitButton("Update")
                      ),
                      tabPanel("Initial prevalence when found",
                               helpText(h4("Set initial prevalence using rule of thumb")),
                               helpText(h5("Currently implemented via scenarios")),
                               tags$hr(),
                               submitButton("Update")
                      ),
                      tabPanel("Scenarios",
                               helpText(h4("Setting scenarios")),
                               tags$hr(),
                               h5("Current parameters:"),
                               h4(tableOutput("currentscenarioTable.scenarios")),
                               tags$hr(),
                               uiOutput("selectScenario"),
                               tags$hr(),
                               h5("Available scenarios:"),
                               h4(tableOutput("scenarioTable.scenarios")),
                               tags$hr(),
                               submitButton("Update")
                      )
           ),
## input and output menu ----
           navbarMenu("Input and output",
                      tabPanel("Manage scenarios",
                               p("Use this tab to save an Excel spreadsheet and subsequently to restore the simulations"),
                               fileInput("uploadPars", "Choose CSV File",
                                           accept = c(
                                             "text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")
                                 ),
                                 tags$hr(),
                                 h5("Available scenarios:"),
                                 h4(tableOutput("scenarioTable.upload")),
                               tags$hr(),
                               h5("Current parameters:"),
                               h4(tableOutput("currentscenarioTable.upload")),
                               tags$hr(),
                                h5("Save current scenario:"),                              
                                downloadButton('downloadPars', 'Download scenario'),
                               tags$hr(),
                               p("Note: Only columns 10 to 24 in the Excel file are important; please note that no checking of parameters is currently performed. The columns represent:"),
                               p("Column 10: Area (default is 100 to represent %)"),
                               p("Column 11: Rate of spread (per year)"),
                               p("Column 12: Initial infected area"),
                               p("Column 13: Value: Timber (once at harvest)"),
                               p("Column 14: Value: Recreation (continuous)"),
                               p("Column 15: Value: Landscape (continuous)"),
                               p("Column 16: Value: Biodiversity (continuous)"),
                               p("Column 17: Value: Carbon (continuous)"),
                               p("Column 18: Control efficacy"),
                               p("Column 19: Control cost"),
                               p("Column 20: Control variability model (0 means no variability, 1 is exponential)"),
                               p("Column 21: Control option (0 all area, 1 only infected area)"),
                               p("Column 22: Number of replicates"),
                               p("Column 23: Discount rate"),
                               p("Column 24: Effort"),
                               tags$hr(),
                               p("Sample spreadsheet (input.csv) can be found",a(href="https://app.box.com/s/t4f08okmj08137r0gn6bhwxdgfuqjrlo","here.")," Please download it to your computer and then use the above menu to upload it to the FC Tool")
                      ),
                      tabPanel("Download results",
                               h5("Choose a dataset (press Update after selection and only then press Download):"),
                               fluidRow(
                                 column(4,selectInput("dataset", NULL, 
                                           choices = c("Healthy", "Infected"))),
                               column(3,submitButton("Update choice"))),
                               downloadButton('downloadData', 'Download data')
                      ),
                      tabPanel("Generate report",
                               fluidPage(
                                 h3("Simulation report"),
                                 tags$hr(),
                                 fluidRow(
                                   helpText(h4("Predictions:")),
                                   column(6,
                                          h5("Proportion of healthy and infected trees 
                                             and total damage over a year:"),
                                          h4(tableOutput("table.report"))
                                          ),
                                   column(6,
                                          h5("Distribution of combined damage:\nred=with control, blue=do nothing"),
                                          plotOutput("distPlot.report",width="275px",height="200px")
                                   )
                                 ),
                                 tags$hr(),
                                 h4("Current list of key parameters:"),
                                 h3(tableOutput("parametersTable.report"))
                              ) # fluid page
           )
           ),
## help menu ----
           tabPanel("Help",
                    fluidPage(
                      withMathJax(helpText(h2("Forest simulator ."))),
                      h4("The description below does not fully describe changes in the model"),
                      hr(),
                      p("The model:"),         
                      hr(),
                      p("$$\\frac{dS}{dt}=-\\beta I\\left(1- \\frac{I}{N}\\right)$$"),
                      p("$$\\frac{dI}{dt}=\\beta I\\left(1- \\frac{I}{N}\\right) - C I$$"),
                      hr(),
                      p("where \\(S\\) represents healthy forest area, \\(I\\) represents infected area, and \\(N=S+I\\)."),
                      p("\\(\\beta\\) is the rate of spread and \\(C\\) represents control efforts."),
                      p("Time is measured in years and the horizon is assumed to be a year. Because of a short period, no discounting is applied."),
                      p("These equations are solved numerically twice, once with and once without control effort."),
                      p("Control effort is a function of investment:"),
                      hr(),
                      p("$$C=c_r e $$"),
                      hr(),
                      p("where \\(c_r\\) is a rate at which control is effective and \\(e\\) is the effort (represented here as the number of man-week)"),
                      p("The total damage as a function of time is given by:"),
                      hr(),
                      p("$$\\Pi(t)=v I(t) + \\int_0^t w_e e \\, dt$$"),
                      p("or"),
                      p("$$\\Pi(t)=v I(t) + \\int_0^t w_e e\\, \\frac{I(t)}{N}\\, dt$$"),
                      hr(),
                      p("where \\(e\\) is the effort (as above) and \\(w_e\\) is the cost of each man-week. \\(v\\) is the value of healthy forest and \\(S(t)\\) is the current number of healthy trees."),
                      p("Additionally, the state is summarised at the end of the year and final share of healthy and infected areas, the final damage and the cumulative cost are all shown."),
                      hr()
                      )
           )
)

## required for running ----
shinyApp(ui = ui, server = server, options = list(launch.browser=T))
