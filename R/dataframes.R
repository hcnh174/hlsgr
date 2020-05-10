#' Load data frame
#'
#' @param filename
#' @param idcol
#' @param header
#' @param stringsAsFactors
#' @param check.names
#' @param comment.char
#' @param sep
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
loadDataFrame <- function(filename, idcol=NULL, header=TRUE, stringsAsFactors=FALSE, check.names=TRUE, comment.char='#', sep='\t', ...)#default.stringsAsFactors())
{
  checkFileExists(filename)
  dataframe <- read.table(filename, header=header, encoding='UTF-8', sep=sep, comment.char=comment.char,
                          stringsAsFactors=stringsAsFactors, check.names=check.names, ...)
  if (!is.null(idcol))
    rownames(dataframe) <- dataframe[[idcol]]
  return(dataframe)
}

######################################################################

#' Split fields on delimiter
#'
#' @param str
#' @param delimiter
#'
#' @return
#' @export
#'
#' @examples
#' splitFields('abc, def, ghi')
splitFields <- function(str, delimiter=',')
{
  if (!is.character(str))
    return(as.array(str))
  if (length(str)==0)
    return(character(0))
  if (length(str)>1)
    return(str)
  if (length(nchar(str))>1)
    return(str)
  fields <- strsplit(str,delimiter)[[1]]
  return(trim(fields))
}

######################################################################

#' Append values to vector
#'
#' @param values1
#' @param values2
#'
#' @return
#' @export
#'
#' @examples
#' appendValues('a,b,c,d,e,f,g','h,i,j,k,l')
#' appendValues(c(),'q')
appendValues <- function(values1,values2)
{
  return(c(splitFields(values1),splitFields(values2)))
}

######################################################################

#' Remove elements from vector
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
removeElements <- function(x, y)
{
  newvalues <- c()
  for (curval in splitFields(x))
  {
    if (!containsElement(y,curval))
      newvalues <- appendValues(newvalues,curval)
  }
  return(newvalues)
}
#removeElements('a,b,c,d,e,f,g','b,f')

######################################################################

#' Exclude columns from dataframe
#'
#' @param data
#' @param cols
#'
#' @return
#' @export
#'
#' @examples
excludeColumns <- function(data, cols)
{
  cols <- splitFields(cols)
  keep <- removeElements(colnames(data),cols)
  return(data[,keep])
}
#data <- excludeColumns(data,'depth')

######################################################################

#' Add rowname as a column
#'
#' @param data
#' @param colname
#'
#' @return
#' @export
#'
#' @examples
addIdColumn <- function(data, colname='ID')
{
  if (!hasColumn(data,colname))
  {
    ids <- as.data.frame(as.character(rownames(data)))
    colnames(ids) <- colname
    data <- cbind(ids,data)
  }
  return(data)
}
#addIdColumn(mirnas)

######################################################################

#' Checks if a dataframe has a column
#'
#' @param data
#' @param col
#'
#' @return
#' @export
#'
#' @examples
hasColumn <- function(data, col)
{
  return(containsElement(colnames(data),col))
}
#hasColumn(top,'ID')

######################################################################

#' Checks if a vector contains a value
#'
#' @param x
#' @param value
#'
#' @return
#' @export
#'
#' @examples
containsElement <- function(x, value)
{
  for (curval in splitFields(x))
  {
    if (curval==value)
      return(TRUE)
  }
  return(FALSE)
}
#containsElement(splitFields('a,b,c,d,e,f,g'),'d')
#containsElement(splitFields('a,b,c,d,e,f,g'),'q')

######################################################################

#' Checks if a column exists in a dataframe
#'
#' @param data
#' @param colname
#'
#' @return
#' @export
#'
#' @examples
colExists <- function(data, colname)
{
  return(!is.na(match(colname, names(data))))
}
#colExists(counts,'ratio6')

######################################################################

#' Throws an error if a column does not exist in a dataframe
#'
#' @param data
#' @param colname
#'
#' @return
#' @export
#'
#' @examples
checkColExists <- function(data, colname)
{
  if (!colExists(data,colname))
    throw('cannot find col with colname: ',colname)
}
#checkColExists(subjects,'hbv_vs_healthy')

######################################################################

#' Write table to file
#'
#' http://r.789695.n4.nabble.com/Appending-strings-at-the-beginning-of-a-text-file-td901370.html
#'
#' @param table
#' @param filename
#' @param verbose
#' @param row.names
#' @param col.names
#' @param eol
#' @param na
#'
#' @return
#' @export
#'
#' @examples
writeTable <- function(table, filename, verbose=TRUE, row.names=FALSE, col.names=TRUE, eol='\n', na='')
{
  #col.names <- c('id',names(table))
  write.table(table, filename, quote=FALSE, row.names=row.names, col.names=col.names, sep='\t', na=na, eol=eol)
  if (row.names==TRUE)# prepend a column name for the row labels
  {
    fConn <- file(filename, 'r+')
    lines <- readLines(fConn)
    lines[1] <- concat('data\t', lines[1])
    writeLines(lines, con=fConn)
    close(fConn)
  }
  if (verbose)
    print(paste('wrote table to file:',filename))
}
