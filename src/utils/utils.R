util = new.env(hash=T)

##NB the append function seems to do wierd things when additional is a dataframe.
util$appendToList <- function(existingList, additional)
{
    existingList[[length(existingList)+1]]=additional
    return(existingList)
}

#intentionally in global scope for ease of typing
fp = file.path
