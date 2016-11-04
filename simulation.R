#set parameters
	N   <- 40   # amount of words
	S   <- 15   # amount of signals
	P   <- .26  # probability of signals
	M   <- 15   # 2^M is max number of languages in a family
	Bsc <- .15    # bursting coefficient for splitting causing change
	Bp  <- 5    # time steps for which bursting persists
#	st  <- 200  # number of time steps
	pl  <- .425   # probability of a lexical change
	pp  <- .425   # probability of a phonological change
	ps  <- .03  # base probability of a split
	pe  <- .004  # probability of extinction

#wrapper function for running simulation given different parameter settings
#X is the number of families created
run <- function(X) {
	for (i in 1:X) {
		steps <- sample(c(50:200), 1, replace = T)
		cat("\nsteps:",steps,"\n")
		fam(topol(M),steps, pl, pp, N, S, P)
	}
}

#define phonological inventory to draw from, distinguishing consonants and vowels
invent <- c("p", "b", "f", "v", "m", "w", "8", "4", "t", "d", "s", "z", "c", "n", "r", "l", "S", "Z", "C", "j", "T", "5", "y", "k", "g", "x", "N", "q", "G", "X", "h", "7", "L", "i", "e", "E", "3", "a", "u", "o")
con <- invent[1:33]
vow <- invent[34:40]

#create an artifical list of N words with:
#a tendency for consonant-vowel structure,
#a probability for different symbols to occur based on the frequency distribution in the world's languages
#word lengths similar to averages across the world's languages
#a number (S) of the words carry some sound-meaning relationships (signals) with
#a certain probability, P (S = 15 and P = .26 are recommended)
ple <- c(.0075, .1077, .2015, .2529, .1817, .1207, .0635, .0343, .0159, .0143) # probability vector for different word lengths
pwo <- c(2.53, 1.95, 0.42, 0.42, 4.61, 2.49, 0.12, 0.002, 4.38, 1.84, 2.45, 0.23, 0.94, 6.04, 2.73, 2.36, 0.92, 0.11, 1.14, 0.21, 0.12, 0.30, 2.27, 5.20, 1.39, 0.79, 0.83, 0.47, 0.03, 0.17, 2.70, 1.57, 0.45, 10.52, 5.96, 0.75, 2.21, 15.42, 6.25, 6.62)
sigs <- c("l", "k", "p", "n", "n", "N", "g", "l", "m", "r", "k", "x", "t", "r", "n", "p", "b", "f", "v", "m", "w", "8", "4", "t", "d", "s", "z", "c", "n", "r", "l", "S", "Z", "C", "j", "T", "5", "y", "k", "g", "x", "N", "q", "G", "X", "h", "7", "L", "i", "e", "E", "3", "a", "u", "o")

wls <- function (N, S, P) {
	words <- c()
	for (i in 1:N) {
		l <- sample(c(1:10), 1, replace = T, prob = ple)
		wch <- c()
			for (j in 1:l) {
				if ( (j %% 2)!=0 ) {
					wch[j] <- sample(invent, 1, prob = pwo*c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.5,0.5,0.5,0.5,0.5,0.5,0.5))
				}
				if ( (j %% 2)==0 ) {
					wch[j] <- sample(invent, 1, prob = pwo*c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2))
				}
			}
		w <- paste(wch, collapse="")
		words[i] <- w
	}
	words <- epen(words)
	words <- epen(words)
	words <- vredux(words)
	words <- vredux(words)
	words <- vadd(words)
	words <- seqred(words)
	for (k in 1:S) {
		toss <- sample( c(1:0), 1, prob=c(P, (1-P) ) )
		if (toss==1) {
			a1 <- unlist(strsplit(words[k],""))
			a2 <- length(a1)
			a3 <- sample(c(1:a2), 1)
			a1[a3] <- sigs[k]
			words[k] <- paste(a1, collapse="")
		}
	}
	return(words)
	rm(words); rm(i); rm(l); rm(wch); rm(j); rm(w); rm(S); rm(k); rm(toss); rm(a1); rm(a2); rm(a3)
}

#epenthesis, inserting a vowel in the middle of a sequence of four consonants
#x is a word list
epen <- function(x) {
	for (i in 1:length(x)) {
		nc <- nchar(x[i])
		if ( nc >= 4 ) {
			for (j in 1:(nc-3)) {
				a1 <- substring(x[i], j, j+3)
				a2 <- unlist(strsplit(a1,""))
				a3 <- a2[!is.na(match(a2,con))]
				if ( length(a3)==4 ) {
					a4 <- unlist(strsplit(x[i],""))
					first <- a4[1:(j+1)]
					second <- a4[(j+2):length(a4)]
					vowel <- sample(vow, 1)
					a5 <- c(first, vowel, second)
					x[i] <- paste(a5, collapse="")
					break
				}
			}			
		}
	}
	return(x)
	rm(x); rm(a1); rm(a2); rm(a3); rm(a4); rm(a5)
}

#add a consonant to the end of a word only consisting of consonants
vadd <- function(x) {
	for (i in 1:length(x)) {
		a1 <- unlist(strsplit(x[i],""))
		a2 <- a1[!is.na(match(a1,con))]
		if ( length(a1)==length(a2) ) {
			vowel <- sample(vow, 1)
			x[i] <- paste(c(a1,vowel),collapse="") 
		}
	}			
	return(x)
	rm(x); rm(vowel); rm(a1); rm(a2)
}

#reduce a sequence of four vowels by deleting the second
vredux <- function(x) {
	for (i in 1:length(x)) {
		nc <- nchar(x[i])
		if ( nc >= 4 ) {
			for (j in 1:(nc-3)) {
				a1 <- substring(x[i], j, j+3)
				a2 <- unlist(strsplit(a1,""))
				a3 <- a2[!is.na(match(a2,vow))]
				if ( length(a3)==4 ) {
					a4 <- unlist(strsplit(x[i],""))
					first <- a4[1:(j)]
					second <- a4[(j+2):length(a4)]
					a5 <- c(first, second)
					x[i] <- paste(a5, collapse="")
					break
				}
			}			
		}
	}
	return(x)
	rm(i); rm(nc); rm(j); rm(a1); rm(a2); rm(a3); rm(a4); rm(a5); rm(x)
}

#reduce a sequence of three identical segments to two
seqred <- function(x) {
	for (i in 1:length(x)) {
		nc <- nchar(x[i])
		if ( nc >= 3 ) {
			for (j in 1:(nc-2)) {
				a1 <- substring(x[i], j, j+2)
				a2 <- unlist(strsplit(a1,""))
				if ( a2[1]==a2[2] ) {
					if ( (a2[2]==a2[3]) ) {
						a3 <- unlist(strsplit(x[i],""))
						a3 <- a3[-j]
						x[i] <- paste(a3, collapse="")
						break
					}
				}
			}			
		}
	}
	return(x)
	rm(x); rm(nc); rm(j); rm(a1); rm(a2); rm(a3); rm(i)
}

#creates a vector of sound change probabilities from Brown et al (2013)
#but C:V correspondences are completely excluded
cor <- read.table(file="correspondences.txt",header=T)
m <- matrix(0, nrow=40, ncol=40)
rownames(m)=invent
colnames(m)=invent
for (i in 1:length(invent)) {
	for (j in 1:length(invent)) {
		c1 <- which(cor[,1]==invent[i])
		c2 <- which(cor[,2]==invent[j])
		if ( length(intersect(c1,c2)) > 0 ) {
			m[i,j] <- cor[intersect(c1,c2),3]/100
			m[j,i] <- cor[intersect(c1,c2),3]/100
		}
	}
}

#do a phonological change affecting a random sound, changing 
#all sounds, a, into another, b, with a probability of 0.5; 
#b is chosen with a probability corresponding to CP
#in Brown et al (2013), but C:V changes are excluded
phon_change <- function(x, N) {
	poslist <- sample(c(1:N), 1)        # choose a random word
	s <- unlist(strsplit(x[poslist],"")) # choose a random sound
	howmany <- length(s)
	posword <- sample(c(1:howmany), 1)
	sound <- s[posword]
	posinvent <- grep(sound, invent)     # find the sound in the inventory
	pt <- m[,posinvent]                  # assign p's for targets of change
	new <- x
	target.no <- sample(c(1:length(invent)), 1, prob = pt)
	target <- invent[target.no]
	for (i in 1:length(new)) {
		split <- unlist(strsplit(new[i],""))
		where <- grep(sound,split)
		if (length(where) > 0) {
			for (j in 1:length(where)) {
				if ( sample(c(0:1),1)==1 ) {
					split[where[j]] <- target
				}
			}
		}
		new[i] <- paste(split, collapse="")
	}
	new <- seqred(new)
	return(new)
}

#lexical change
lex_change <- function(x) { # x is a wordlist
	old.word <- sample(c(1:length(x)),1)
	new.word <- wls(1, 1, .080859375)
	x[old.word] <- new.word
	return(x)
}

#make table guiding tree-growing in binarily splitting family
topol <- function(M) {
	bigsis <- c(1)
	littlesis <- c()
	for (i in 1:M) {
		bigsis <- c(bigsis, 1:2^(i-1))
	}
	littlesis <- c(2:length(bigsis))
	bigsis <- bigsis[-1]
	l <- length(bigsis)
	guide <- as.matrix(cbind(bigsis, littlesis, c(rep(NA,l))))
	return(guide)
}

#prune the topology scheme in case of extinction or non-splitting
#a: dying lineage
#b: topology scheme
#at: place in the scheme where extinction happens,
#this can either be just before a split (left in topology scheme)
#or just after (right); if before a split the program
#backtracks to the birth and erases that and it turns offspring
#into offspring of the dying lineage's mother; if after a split
#birth and all offspring are simply removed
death <- function(a,b,at) {
	pos <- ""
	if ( a==b[at,1] ) {
		pos <- "left"
	}
	if ( a==b[at,2] ) { # if da dies equivalent to extinction of mo's descendants
		a <- b[at,1]
	}
	e <- c(a, b[at,2])
	b <- b[-at,,drop=FALSE]
	if ( at < length(b[,1]) ) {
		i <- at - 1
		kill <- c()
		k <- 0
		for ( i in at:length(b[,1]) ) { # kill descendants
			if ( length(intersect(e,b[i,1])) > 0 ) {
				e <- c(e, b[i,2])
				k <- k + 1
				kill <- c(kill,i)
			}
		}
		if ( length(kill) > 0 ) {
			b <- b[-kill,,drop=FALSE]
		}
	}
	if ( pos=="left" & length(b[,1]) > 0 ) { # left position, i.e., extinction
		if ( a==1 ) {
			h <- min(match(a,b[,1])) # find birthplace
			g <- b[h,2] # identify mother
		}
		if (a!=1 ) {
			h <- min(match(a,b[,2])) # find birthplace
			g <- b[h,1] # identify mother
		}
		if (a!=1 & is.na(h)) {
			h <- min(match(a,b[,1])) # find birthplace
			g <- b[h,2] # identify mother
		}
		if ( length(b[,1]) > 0 ) {
				b <- b[-h,,drop=FALSE] # removes birth
		}
		if ( !is.na(match(a,b[,1])) ) { # updates previous births to be from mother
			j <- which(a==b[,1])
			for (k in 1:length(j)) {
				b[j[k],1] <- g
			}
		}
	}
	return(b)
}

#create word lists for the stages in a linguistic lineage
stepslex <- function(st, pl, pp, proto, N, Bsc, Bp) {
	stages <- list()
	changes <- list()
	stages[[1]] <- proto
	cf <- read.table(file="currentfam.txt")[1,1]
	cs <- read.table(file="currentstep.txt")[1,1]
	for (i in 1:(st-1)) {
		changes[[i+1]] <- c(0,0,0,0,cf,i+cs-1) # bpc, blc, npc, nlc, lineage, step
		if ( i <= Bp ) {
			lch <- sample(c(1,0),1,prob=c(pl+Bsc,(1-pl-Bsc)))
			pch <- sample(c(1,0),1,prob=c(pp+Bsc,(1-pp-Bsc)))
			if (pch==1) {
				stages[[i+1]] <- phon_change(stages[[i]],N) 
				changes[[i+1]][1] <- 1
			}
			if (pch==0) {
				stages[[i+1]] <- stages[[i]]
			}
			if (lch==1) { 
				stages[[i+1]] <- lex_change(stages[[i+1]]) 
				changes[[i+1]][2] <- 1
			}
			if (lch==0) {
				stages[[i+1]] <- stages[[i+1]]
			}
		}
		if ( i > Bp ) {
			lch <- sample(c(1,0),1,prob=c(pl,(1-pl)))
			pch <- sample(c(1,0),1,prob=c(pp,(1-pp)))
			if (pch==1) { 
				stages[[i+1]] <- phon_change(stages[[i]],N) 
				changes[[i+1]][3] <- 1
			}
			if (pch==0) {
				stages[[i+1]] <- stages[[i]]
			}
			if (lch==1) { 
				stages[[i+1]] <- lex_change(stages[[i+1]]) 
				changes[[i+1]][4] <- 1
			}
			if (lch==0) {
				stages[[i+1]] <- stages[[i+1]]
			}
		}
	}
	for (j in 1:length(changes)) {
		cat(unlist(changes[[j]]), "\n", file="changes.cha", sep="\t", append=T)
	}
	return(stages)
}

#takes the output of stepslex, introduces a split with a certain probability
#and if there is no split it returns 0, and descendants should then be removed
#ps is base probability of a split, pe is extinction probability
stepssplit <- function(x,ps,pe) {
	for (i in 1:(length(x)-1)) {
		toss <- sample(c(-1,1,0),1,prob=c(pe, ps, (1-pe-ps)))
		if ( toss==1 ) {
		return(i)
		break
		}
		if ( toss==-1 ) {
		return(-1)
		break
		}
	}
	return(0)
}

#output terminal taxa in ASJP format in a *.dat file and parameter settings in a *.par file
TO <- function(F, st, N, fn) {
	term <- list()
	names <- c()
	c2 <- 0
	for (i in 1:length(F)) {
		if ( length(F[[i]]) > 1 ) {
			c2 <- c2 + 1
			term[[c2]] <- unlist(F[[i]][st])
			names <- c(names,i)
		}
	}
	
	gl <- c("1 I\t", "2 you\t", "3 we\t", "11 one\t", "12 two\t", "18 person\t", "19 fish\t", "21 dog\t", "22 louse\t", "23 tree\t", "25 leaf\t", "28 skin\t", "30 blood\t", "31 bone\t", "34 horn\t", "39 ear\t", "40 eye\t", "41 nose\t", "43 tooth\t", "44 tongue\t", "47 knee\t", "48 hand\t", "51 breast\t", "53 liver\t", "54 drink\t", "57 see\t", "58 hear\t", "61 die\t", "66 come\t", "72 sun\t", "74 star\t", "75 water\t", "77 stone\t", "82 fire\t", "85 path\t", "86 mountain\t", "92 night \t", "95 full\t", "96 new\t", "100 name\t")
	fnd <- paste(fn,".dat",sep="")
	pa <- readLines("preamble.txt")
	cat(pa,file=fnd,sep="\n")
	for (j in 1:length(names)) {
	cat(paste("L",names[j],"{Fam.GENUS|Ethnologue@Glottologue}\n 3  100.00  100.00       10000         iso\n",sep=""),file=fnd,append=T)
		for (k in 1:N) {
			cat(paste(gl[k],unlist(term[[j]])[k],"\n",sep=""),file=fnd,append=T)
		}
	}
	cat("     \n",file=fnd,append=T)
	fnp <- paste(fn,".par",sep="")
	cat("N =",N,"amount of words\n",file=fnp)
	cat("S =",S,"amount of sound symbolic signals\n",file=fnp, append=T)
	cat("P =",P,"probability of signals\n",file=fnp, append=T)
	cat("M =",M,"2^M is max number of languages in a family\n",file=fnp, append=T)
	cat("Bsc =",Bsc,"bursting coefficient split -> change\n",file=fnp, append=T)
	cat("Bp =",Bp,"time steps during which bursting persists\n",file=fnp, append=T)
	cat("st =",st,"number of time steps\n",file=fnp, append=T)
	cat("pl =",pl,"probability of a lexical change\n",file=fnp, append=T)
	cat("pp =",pp,"probability of a phonological change\n",file=fnp, append=T)
	cat("ps =",ps,"base probability of a split\n",file=fnp, append=T)
	cat("pe =",pe,"probability of extinction\n",file=fnp, append=T)
}

#select the phonological and lexical changes pertaining to just those lineages
#that survive in the final tree and output them to changes.cha, which gets a unique
#name in the fam function
housekeeping.changes <- function(fnt) {
	x <- read.table(file="changes.cha", header=T)
	delete <- c()
	for (i in length(x[,1]):1) {
		m1 <- which(x[1:(i-1),5]==x[i,5])
		m2 <- which(x[1:(i-1),6]==x[i,6])
		m3 <- intersect(m1,m2)
		if ( length(m3) > 0 ) {
			delete <- c(delete, m3)
		}
	}
	x <- x[-delete,]
	y <- read.table(file=fnt, header=FALSE)
	lt <- union(y[,1],y[,2])
	lc <- as.vector(as.numeric(x[,5]))
	d <- setdiff(lc, lt)
	if ( length(d) > 0 ) {
		inds <- which(lc %in% d)
		x <- x[-inds,]
	}
	rs <- min(y[,3])
	br <- which(as.vector(as.numeric(x[,6])) <= rs)
	if ( length(br > 0) ) {
		x <- x[-br,]
	}
	sum1 <- sum(as.numeric(x[,1]))
	sum2 <- sum(as.numeric(x[,2]))
	sum3 <- sum(as.numeric(x[,3]))
	sum4 <- sum(as.numeric(x[,4]))
	perc1 <- round(100*sum1/(sum1+sum3),2)
	perc2 <- round(100*sum2/(sum2+sum4),2)
	write.table(x[,1:6], file="changes.cha", sep="\t", quote=FALSE, row.names=FALSE)
	cat(sum1, "\t", sum2, "\t", sum3, "\t", sum4, "\n", file="changes.cha", append=T)
	cat(perc1, "\t", perc2, "\n", file="changes.cha", append=T)
	return(NULL)
}

#core function for creating a family following the topology scheme, t. The
#topology is put in a *.top file
fam <- function(t, st, pl, pp, N, S, P) {
	cat("\nCREATE NEW FAMILY\n\n")
	cat("\nCREATE NEW FAMILY\n\n",file="log.txt")
	cat("PE_ph\tPE_lx\tNm_ph\tNm_lx\tlineage\tstep\n", file="changes.cha")
	tl <- 0
	F <- list()
	bd <- list()
	proto <- wls(N, S, P)
	write.table(1, file="currentfam.txt")
	write.table(1, file="currentstep.txt")
	F[[1]] <- stepslex(st, pl, pp, proto, N, Bsc, Bp)
	bd[[1]] <- 1
	firststep <- 1
	while ( length(t[,1]) > tl & length(t[,1]) > 0 ) {
		tl <- tl + 1
		Mo <- as.vector(t[tl,1])
		Da <- as.vector(t[tl,2])
		cat("Mo:",Mo,"\t")
		cat("Birth Mo at",bd[[Mo]],"\n")
		cat("Mo:",Mo,"\t",file="log.txt", append=T)
		cat("Birth Mo at",bd[[Mo]],"\n", file="log.txt", append=T)
		sout <- stepssplit( F[[Mo]] [bd[[Mo]]:st], ps, pe)
		s <- sout + bd[[Mo]] - 1
		if ( sout==-1 & firststep==1 ) { 
			t <- topol(1)[-1,,drop=FALSE]
			}
		if ( sout==-1 & firststep!=1 ) { 
			cat("extinction of lineage",Mo,"\n\n")
			cat("extinction of lineage",Mo,"\n\n",file="log.txt", append=T)
			t <- death(Mo,t,tl)
			tl <- tl - 2
			F[[Mo]] <- NA
		}
		if ( sout==0 ) {
			cat("nonsplitting of lineage",Mo,"\n\n")
			cat("nonsplitting of lineage",Mo,"\n\n",file="log.txt", append=T)
			t <- death(Da,t,tl)
			tl <- tl - 1
		}
		if ( s > (st-2) ) {
			cat("attempt to split lineage",Mo,"at tip\n\n")
			cat("attempt to split lineage",Mo,"at tip\n\n",file="log.txt", append=T)
			t <- death(Da,t,tl)
			tl <- tl - 1 
		}
		if ( (sout > 0) & (s < st-1) ) {
			F[[Da]] <- F[[Mo]]
			write.table(Da, file="currentfam.txt")
			write.table(s+1, file="currentstep.txt")
			nDa <- list()
			nDa <- stepslex((st-s+1), pl, pp, unlist(F[[Mo]][s]), N, Bsc, Bp)
			nDa <- nDa[2:length(nDa)]
			F[[Da]][(s+1):st] <- nDa
			write.table(Mo, file="currentfam.txt")
			write.table(s+1, file="currentstep.txt")
			nMo <- list()
			nMo <- stepslex((st-s+1), pl, pp, unlist(F[[Da]][s]), N, Bsc, Bp)
			nMo <- nMo[2:length(nMo)]
			F[[Mo]][(s+1):st] <- nMo
			bd[[Mo]] <- s; bd[[Da]] <- s
			cat("Da:",Da,"\t")
			cat("Birth Da at",bd[[Da]],"\n\n")
			cat("Da:",Da,"\t", file="log.txt", append=T)
			cat("Birth Da at",bd[[Da]],"\n\n", file="log.txt", append=T)
			t[tl,3] <- s
		}
		firststep <- 0
	}
	c0 <- c()
	c1 <- 0
	for (i in 1:length(F)) {
		if ( length(F[[i]]) > 1 ) { # NA lineages are 1 long, null lineages 0 long
			c1 <- c1 + 1
			c0[c1] <- i
		}
	}
	if ( length(t[,1]) > 0 ) {
		firstfull <- min(c0)
		fn <- paste("F_",c1,"_",unlist(F[[firstfull]][st])[1],sep="")
		fnt <- paste(fn,".top",sep="")
		TO(F, st, N, fn)
		for (y in 1:length(t[,1])) {
			cat(t[y,1],"\t",t[y,2],"\t",t[y,3],"\n", file=fnt, append=T)
		}
		housekeeping.changes(fnt)
		file.rename("changes.cha", paste("F_",c1,"_",unlist(F[[firstfull]][st])[1],".cha",sep=""))
	}
	if ( length(t[,1])==0 ) {
		cat("Entire family gone extinct\n\n\n")
		cat("Entire family gone extinct\n\n\n", file="log.txt", append=T)
 	}
	file.remove("currentfam.txt")
	file.remove("currentstep.txt")
	file.remove("log.txt")
	return(NULL)
}

