#This thing requires product 1 and product 2 before we've run the 90% AM+EM cutoff.
tree.table.sub <- tree.table[tree.table$PLT_CN %in% plot.table[plot.table$relEM.AM < 0.9,]$PLT_CN,]
myc.ref <- read.csv('required_products_utilities/mycorrhizal_SPCD_data.csv')
tree.table.sub <- merge(tree.table.sub, myc.ref[,c('SPCD','GENUS','SPECIES')], all.x = T)
weird.trees <- data.frame(table(tree.table.sub$SPCD))
colnames(weird.trees) <- c('SPCD','stems')
weird.trees <- merge(weird.trees,myc.ref, all.x = T)
weird.trees <- weird.trees[order(weird.trees$stems, decreasing = T),]

#which trees in the <90% Am-EM plots are not AM-EM trees, and therefore could be driving this?
#3/4 of the non AM-EM trees are aspen.
test <- weird.trees[!(weird.trees$MYCO_ASSO %in% c('ECM','AM')),]
sum(test[test$GENUS == 'Populus',]$stems, na.rm=T) / sum(test$ste) #73% of stems are populus.

#how many of these plots are driven by populus greater than 10% abundance?
#67% of the non AM-EM stands are driven by populus. Populus is >10% of the trees in the plot for 4037/6025 stands.
#14% of the non-AM-EM stands are driven by trees with unkown mycorrhizal status.
tree.table.sub$populus <- ifelse(tree.table.sub$GENUS == 'Populus',1,0)
tree.table.sub$unk     <- ifelse(tree.table.sub$MYCO_ASSO == 'UNK', 1, 0)
tree.table.sub$basal.unk <- tree.table.sub$BASAL * tree.table.sub$unk
tree.table.sub$basal.pop <- tree.table.sub$BASAL * tree.table.sub$populus
basal <- aggregate(BASAL ~ PLT_CN, data = tree.table.sub, FUN = sum)
basal.pop <- aggregate(basal.pop ~ PLT_CN, data = tree.table.sub, FUN = sum)
basal.unk <- aggregate(basal.unk ~ PLT_CN, data = tree.table.sub, FUN = sum)
test <- merge(basal, basal.pop)
test <- merge(test , basal.unk)
test$rel.pop <- test$basal.pop / test$BASAL
test$rel.unk <- test$basal.unk / test$BASAL
nrow(test[test$rel.pop > 0.10,])
nrow(test[test$rel.unk > 0.10 & test$rel.pop <= 0.10,])

#Checkout the non-populus plots.
tree.table.subpop <- tree.table.sub[tree.table.sub$PLT_CN %in% test[test$rel.pop < 0.1,]$PLT_CN,]
weird.trees <- data.frame(table(tree.table.subpop$SPCD))
colnames(weird.trees) <- c('SPCD','stems')
weird.trees <- merge(weird.trees,myc.ref, all.x = T)
weird.trees <- weird.trees[!(weird.trees$MYCO_ASSO %in% c('ECM','AM')),]
weird.trees <- weird.trees[order(weird.trees$stems, decreasing = T),]
