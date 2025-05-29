# Given a cmap data structure with some number of NA values, naively (2-way) predict all missing values as a function of all non-missing values
predict_missing_naive = function(dat){
	ncells = length(dat)
	cells = names(dat)
	ndrugs = dim(dat[[1]])[[1]]
	drugs = rownames(dat[[1]])

	new_dat = dat

	isNa = matrix(FALSE,nrow=ncells,ncol=ndrugs,dimnames=list(cells,drugs))

	for(i in cells){
		for(j in drugs){
			if(is.na(dat[[i]][j,1])){
				isNa[i,j] = TRUE
			}
		}
	}

	for(i in cells){
		for(j in drugs){
			if(isNa[i,j]){
				otherCells = setdiff(cells,c(i))
				otherDrugs = setdiff(drugs,c(j))
                otherDrugsWithData = intersect(otherDrugs,names(which(!isNa[i,otherDrugs])))
                otherCellsWithData = intersect(otherCells,names(which(!isNa[otherCells,j])))
                M = matrix(0,nrow=0,ncol=dim(dat[[1]])[[2]])
                for(k in otherCellsWithData){
					M = rbind(M,dat[[k]][j,])
                }
                otherCellVector = colMeans(M)
                otherDrugVector = colMeans(dat[[i]][otherDrugsWithData,])
                new_dat[[i]][j,] = colMeans(rbind(otherCellVector,otherDrugVector))
			}
		}
	}

	return(new_dat)
}
