## build functions for model formulae
## generic glm formula construction function; inputs as:
##   resp =>col 1 in dataframe such that r.col=1, 
##   preds=>col 2 thru ncol in dataframe such that p.col=2
##	 spreds=>smoothed terms
##   NOTE: vars as factors must be coerced PRIOR to formula construction
## example call: mod.form(dat1,1,2)


# with quadratic terms
mod.form.quad=function(dat,r.col,p.col) {
  ## generic glm formula construction function; inputs as:
  ##   resp =>col 1 in dataframe such that r.col=1, 
  ##   preds=>col 2 thru ncol in dataframe such that p.col=2
  ##   NOTE: vars as factors must be coerced PRIOR to formula construction
  ## example call: mod.form(dat1,1,2)
  n.col=ncol(dat)                    # No. columns in dataframe
  resp=colnames(dat[r.col])          # assign resp column name
  resp=paste("as.factor(",colnames(dat[r.col]),")",sep="") # assign resp column name
  pred=colnames(dat[c(p.col:n.col)]) # assign preds column names
  pred2=paste("I(",pred,"^2)",sep="")
  mod.formula=as.formula(paste(resp,"~",paste(pred,pred2,sep="+",collapse="+"))) # build formula
  }

# linear terms only
mod.form.lin=function(dat,r.col,p.col) {
  ## generic glm formula construction function; inputs as:
  ##   resp =>col 1 in dataframe such that r.col=1, 
  ##   preds=>col 2 thru ncol in dataframe such that p.col=2
  ##   NOTE: vars as factors must be coerced PRIOR to formula construction
  ## example call: mod.form(dat1,1,2)
  n.col=ncol(dat)                    # No. columns in dataframe
  resp=colnames(dat[r.col])          # assign resp column name
  resp=paste("as.factor(",colnames(dat[r.col]),")",sep="") # assign resp column name
  pred=colnames(dat[c(p.col:n.col)]) # assign preds column names
  mod.formula=as.formula(paste(resp,"~",paste(pred,collapse="+"))) # build formula
  }
  
##default smoothers (df=4)
mod.form.4=function(dat,r.col,p.col) {
  n.col=ncol(dat) # No. columns in dataframe
  resp=colnames(dat[r.col]) # assign resp column name
  resp=paste("as.factor(",colnames(dat[r.col]),")",sep="") # assign resp column name
  pred=colnames(dat[c(p.col:n.col)]) # assign preds column names
  spred = paste("s(", pred, ")", sep="")
  mod.formula=as.formula(paste(resp,"~",paste(spred,collapse="+"))) # build formula
  }

## smoother df=3
mod.form.3=function(dat,r.col,p.col) {
  n.col=ncol(dat) # No. columns in dataframe
  resp=colnames(dat[r.col]) # assign resp column name
  resp=paste("as.factor(",colnames(dat[r.col]),")",sep="") # assign resp column name
  pred=colnames(dat[c(p.col:n.col)]) # assign preds column names
  spred = paste("s(", pred, ",3)", sep="")
  mod.formula=as.formula(paste(resp,"~",paste(spred,collapse="+"))) # build formula
  }

##smoother df=2
mod.form.2=function(dat,r.col,p.col) {
  n.col=ncol(dat) # No. columns in dataframe
  resp=colnames(dat[r.col]) # assign resp column name
  resp=paste("as.factor(",colnames(dat[r.col]),")",sep="") # assign resp column name
  pred=colnames(dat[c(p.col:n.col)]) # assign preds column names
  spred = paste("s(", pred, ",2)", sep="")
  mod.formula=as.formula(paste(resp,"~",paste(spred,collapse="+"))) # build formula
  }



