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
    while (c1_i<=min(e1,e2) && c2_i<=min(e1,e2))
        #find the index of the first character matches seq1[c1_i] in seq2 
        if (occursin(seq1[c1_i],seq2[c2_i:e2])&& c1_i<=min(e1,e2) && c2_i<=min(e1,e2))
            mtch2=findfirst(seq1[c1_i],seq2[c2_i:e2])
            s_s+=1
            c2_i=c2_i+mtch2
            c1_i=c1_i+1
        else
            break
        end
            #for seq1 start from the next position c1_i+1:e1 search seq[c2_i]
            #c2_i is the next position of the matched position 
        if (occursin(seq2[c2_i],seq1[c1_i:e1]) && c1_i<=min(e1,e2) && c2_i<=min(e1,e2))
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


DNA=["CATCTCGGGGCCAAGGTGTC","CGTACAGGCGTTCAGGTGTG","GAATTAGCTATATCGCCCGG","GTGCGCTCCGGCCGCGTCGC"]
d=zeros(Int8,4)
for i=1:length(DNA)
    d[i]=rand(1:length(DNA[i]))
end
d    
sub1=[DNA[1][1:d[1]],DNA[2][1:d[2]],DNA[3][1:d[3]],DNA[4][1:d[4]]]
sub2=[DNA[1][(d[1]+1):length(DNA[1])],DNA[2][(d[2]+1):length(DNA[2])],DNA[3][(d[3]+1):length(DNA[3])],DNA[4][(d[4]+1):length(DNA[4])]]
sub=[]
append!(sub,[sub1])
append!(sub,[sub2])
seq1=sub[1][1]
seq2=sub[1][2]
sc=Pairwise_Score(seq1,seq2)
sc