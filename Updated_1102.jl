using Random

#longest string in DNA subseqs
function Longest_Seq(DNA)
    #length of DNA seqs
    len=zeros(Int16,length(DNA))
    f=0
    for i=1:length(DNA)
        len[i]=length(DNA[i])
        f=findmax(len)[2]
    end
    return f
end

#PairSum--approximation of alignment score
#input two sequences
#output sum of score 
#from position 1 to length of seq1/2
function Pairwise_Score(seq1,seq2)
    e1=length(seq1)
    e2=length(seq2)
    c1_i=1
    c2_i=1
    s_s=0
    while (c1_i<=e1 && c2_i<=e2)
        #find the index of the first character matches seq1[c1_i] in seq2 
        if (occursin(seq1[c1_i],seq2[c2_i:e2]))
            mtch2=findfirst(seq1[c1_i],seq2[c2_i:e2])
            s_s+=1
            c2_i=c2_i+mtch2
            c1_i=c1_i+1
        else
            break
        end
            #for seq1 start from the next position c1_i+1:e1 search seq[c2_i]
            #c2_i is the next position of the matched position 
        if ( c1_i<=e1 && c2_i<=e2 && occursin(seq2[c2_i],seq1[c1_i:e1]))
            mtch1=findfirst(seq2[c2_i],seq1[c1_i:e1])
            s_s+=1
            c1_i=c1_i+mtch1
            c2_i=c2_i+1
        else
            break
        end
    end
    return s_s
end


#cut subsequence for a d vector
function Generate_Subseq(DNA,d)
    sub=[]
    sub1=[]
    sub2=[]
    for t=1:length(DNA)
        append!(sub1,[DNA[t][1:d[t]]])
        append!(sub2,[DNA[t][(d[t]+1):length(DNA[t])]])
    end
    append!(sub,[sub1])
    append!(sub,[sub2])
    return sub
end

#score of subsequences generated from 1 d_cut position vector
function Cut_Score(DNA,d)
    scr_d=0
    sub=Generate_Subseq(DNA,d)
    for s=1:length(sub)
        for i=1:(length(sub[s])-1)
            for j=(i+1):length(sub[s])
                seq1=sub[s][i]
                seq2=sub[s][j]
                scr_d+=Pairwise_Score(seq1,seq2)
                scr_d+=Pairwise_Score(seq2,seq1)
                scr_d+=Pairwise_Score(reverse(seq1),reverse(seq2))
                scr_d+=Pairwise_Score(reverse(seq2),reverse(seq1))
            end
        end
    end
    return scr_d
end               

#Randomly generate a population of k=100, each represent a possible solution of a cutoff
#DNA is the DNA sequences for multiple sequence alignment
#K is the number of possible solutions 
function Generate_kPop(DNA,K)
    pop=[]
    for k=1:K
        #1 possible solution
        d=zeros(Int16,length(DNA))

        for i=1:length(DNA)
            d[i]=rand(1:length(DNA[i]))
        end
        append!(pop,[d])
    end
    return pop
end    

#all k possible subseqs(sub1+sub2)
function K_Subseqs(DNA,k_pop)
    subseq=Matrix(undef,1,length(k_pop))
    for k in 1:length(k_pop)
        d=k_pop[k]
        subseq[k]=Generate_Subseq(DNA,d)
    end
    return subseq
end

#the next population: new person in a k_pop
function Generate_new_p(DNA,k_pop,idx_best)
    #the longest seuqence f
    f=Longest_Seq(DNA)
    new=zeros(Int16,length(k_pop[idx_best]))
    for i in 1:length(new)
        new[i]=k_pop[idx_best][i]-k_pop[idx_best][f]+floor(1/2*length(DNA[f]))
    end  
    return new
end

#score of all k possible subseqs(sub1+sub2)
function K_Scores(DNA,k_pop)
    len=length(k_pop)
    d_score=zeros(Int16,len)
    k_subs=K_Subseqs(DNA,k_pop)
    for k in 1:length(k_subs)
        d=k_pop[k]
        d_score[k]=Cut_Score(DNA,d)
    end
    return d_score 
end

function Cut_Best(DNA, k, Cycle, T)
    #the longest seuqence f
    f=Longest_Seq(DNA)
    # the founders in a population with k individual possible solutions 
    k_pop=Generate_kPop(DNA,k)
    k_scr=zeros(Int16,k)
    idx_best=0
    scr_best=0
    idx_wrst=0
    scr_wrst=0
    
    for gen=1:Cycle
        k_scr=K_Scores(DNA,k_pop)
        
        #find the current best score/current worst score
        idx_best=findmax(k_scr)[2]
        scr_best=findmax(k_scr)[1]
        idx_wrst=findmin(k_scr)[2]
        scr_wrst=findmin(k_scr)[1]

        #generate a new cutoff d-new, 
        d_new=Generate_new_p(DNA,k_pop,idx_best)
        scr_new=Cut_Score(DNA,d_new)
        #compare d_new with the worst socre
        #update the datasets k_pop,k_score
        if scr_new>scr_wrst
            k_pop[idx_wrst]=d_new
            k_scr[idx_wrst]=scr_new
        end
        #find the updated best/worst score
        idx_best=findmax(k_scr)[2]
        scr_best=findmax(k_scr)[1]
        idx_wrst=findmin(k_scr)[2]
        scr_wrst=findmin(k_scr)[1]

        # generate new cutoff by crossover and mutation
        for t in 1:T
            d_new=zeros(Int16,length(DNA))
            #chose two individual for X&M
            mt=rand((1:length(k_pop)),2)
            i=mt[1]
            j=mt[2]

            #mutation
            #which DNA mutate, cut position  
            c=rand(1:length(DNA))
            #mutate to a new cut position
            ran=rand(1:length(DNA[c]))

            #mutation and Xrsover in seq1
            d_new=[k_pop[i][1:c];k_pop[j][(c+1):length(k_pop[j])]]
            #a mutation happens here
            d_new[c]=ran

            scr_new=Cut_Score(DNA,d_new)
            #update the datasets k_pop
            if scr_new>scr_wrst
                k_pop[idx_wrst]=d_new
                k_scr[idx_wrst]=scr_new
            end
            #update best/ worst score
            idx_best=findmax(k_scr)[2]
            scr_best=findmax(k_scr)[1]
            idx_wrst=findmin(k_scr)[2]
            scr_wrst=findmin(k_scr)[1]

            #X&M in mc2
            d_new=[k_pop[j][1:c];k_pop[i][(c+1):length(k_pop[i])]]
            scr_new=Cut_Score(DNA,d_new)
            #update the datasets k_pop
            if scr_new>scr_wrst
                k_pop[idx_wrst]=d_new
                k_scr[idx_wrst]=scr_new
            end
            #update best/ worst score
            idx_best=findmax(k_scr)[2]
            scr_best=findmax(k_scr)[1]
            idx_wrst=findmin(k_scr)[2]
            scr_wrst=findmin(k_scr)[1]
        end
        #print(k_pop[idx_best])
    end
    return idx_best,k_pop[idx_best],k_scr[idx_best]
end

#read in test sequences
using FastaIO
sub_DNA=[] 
FastaReader("test1_seq.fasta") do fr
    for (desc, seq) in fr
        append!(sub_DNA,["$seq"])
    end
end
sub_DNA=[sub_DNA]

#cut test sequence
cut=Cut_Best(sub_DNA,40, 100, 100)
sub=Generate_Subseq(sub_DNA,cut[2])

function Cut_Subs(sub,k, Cycle, T)
    new_sub=[]
    for i=1:length(sub)
        f=Longest_Seq(sub[i])
        if length(sub[i][f])<=10
            append!(new_sub,sub[i])
        else
            cut=Cut_Best(sub[i],k, Cycle, T)
            append!(new_sub,Generate_Subseq(sub[i],cut[2]))
        end
    end
    return new_sub
end


sub=Cut_Subs(sub_DNA,40, 100, 100)
sub=Cut_Subs(sub,40, 100, 100)


8-element Vector{Any}:
 Any["AACATTGATGCCACATAGTGGTCAGAACAATAGTCGTGTGA", "AACATCGATGCCACATAGTGGTCGAAATAATGGTCTGAGG", "AACCGTGAGCTTGCAGACAGAACACTCCTACTCTGGTGG", "TTACTTCATTCACACCAGAACACAGCCTTTCTTAGATGTCTAGGA", "TGGTGGTAGTGACTCTCTAGACTATAACAAGTGCT", "TTCTGTTCGTGAAGCAAACCGCATTATATGACACTACTTATCAGG", "TGCATATACCAAGTTACCCCGCATTATACGCCTCTT", "ATACTATAGCGAGCTCATGATATAAAGTTGTTAACAGGTTAT"]
 Any["AAGGGTCATCTCTC", "GAACGCAAGTCTCC", "GGCGGACTCTCAAGGT", "CATTCTCTCTGTTCTAATCAC", "TCTTTCTAAAGACAA", "AACAAGATCCTCACC", "CTTCCCAAGGACCAAATCCTAACCATTCT", "GGTATAGCAATTCCCCTACCATTCAC"]
 Any["GACTGCTTTGCACTGAATGAG", "CGACTGCATTGCGCACTGAATAAGAG", "TATGTACTGCCTCAGACGGAGACCAG", "AGCTGGGACAGTCACCTGCTTAG", "AAGAGCCACCTGGCCCAGAAGGGCGGGAGTAGAGCATCTCTTT", "AGCCTCCTATAAATGGGGTACGGCATAGATTGCGG", "GTTGTAGATCGGGGACGGAGTTTCGGTATCT", "ACGTAACCACTCGCAACCCCGCGTGGT"]
 Any["AGAGGAACCTGACAGTACTTCTCTGAGTGGACGACGGGCGTT", "AGGCAATAGTCGCAGTATTCTCCCGTG", "CGCGGAGCCGTCACCGTGACCTGCGGTCGT", "GCCGCTCACGCCGAACCACATTACAATGTACGAC", "ATATTAGGATGGCAACGGACGCGGACGGTTTAGAGGAAGAA", "GAGAGGGGCGGTAATGGACCGGTCCTGGGTCTAGAAGCGGCCCG", "TTAGGATCACGAATGACTTGTCTGTTCCCTTGG", "GTGGTCAGGTTTGATGCGATTAGACGAACGCC"]
 Any["CATCGAACCCACGCCAGTATTCGCGAAC", "GACAACGGGCTTTCATATC", "TCATGGAGCATGCGC", "TCAAGCCAGACTACGGGCATAGCGT", "TGACTTATTGACTTGT", "TCGCACTCGAC", "TTGGTCCCCCCGGGTG", "CGGAGCGAAGCATTTAACTG"]
 Any["ACCGACTCTGCCGATTGCTTAGTCCTTATTAATGCGTTACGA", "ACACTATTCGCGCACTCCGGCTCTGCAGACTGCTTACTCCTTACTAACG", "ACGTAATCGCTGTCCCCGATTCTTCTGTGTGCTTAGGCTGTCGGGTTCGGTTATCACG", "CTCGGCCAGTATTTGATACCAACGAGTGTCATTCTCCCGGCGCCCAGATCGAAA", "TATTCTAGCTACTGAGCAGACGCCTTTATAGTCAACCCTCACGTAGATG", "AAATCTACTTCGCTTATAGGTTCACATATCCGTGCGGACGGTGGAA", "CTTGTCTATAGACTCGCAGGTCTTATAGGCTCCCAAAGGGGAG", "CAGCTTCTGCATACTTGAGCACTAACTGCGTTCGGTTCATCAGTA"]
 Any["CGTGTCCAAAGAGAGCAATTTAATAATCGCGTTGCATCTAGGGT", "CGCTACAACAGGTCCAAAGAAGACAATACTACCCTCC", "TTCCAAATGACAGGAGATAAAGGATAGCAGTACATTATGAGGT", "GCTACAGCGGCCAGCCTTTTATGTAACGCTGCC", "TGTTCTGCAAGCTGAGACGGTGCGTC", "AAAGTGTTCTGGAAGTATGCGGGTAGGATATCAGCC", "CGCGAAGTTTTTATATTCAAGAGGAAGTAGGCGGGAATCTC", "TGTGAGTGCGAGCGTGTGCTCCGCCATAGGGGGGGGGAGTTAAAAGTCGCCC"]
 Any["GTTATCGATCTCGGACAGGGT", "AAGCGTGTATGCGGTCTAGGGCAGAGA", "AGGTGATTGAACATCGGGCGGCT", "ATAGATTGAGCGGCAACCTGAT", "ATATTGATCCCCCACGTTGTGGTGTC", "ATATAACGTTACCTG", "ATAGTACAAAAACGCTTGTATT", "GTTATCATACTATGCCAGGTGAT"]