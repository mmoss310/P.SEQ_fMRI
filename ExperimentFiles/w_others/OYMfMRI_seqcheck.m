clc

%SEQUENCE SETTING
SEQ_ST.BASE(1,:)=[1,2,3];
SEQ_ST.BASE(2,:)=[1,3,2];
SEQ_ST.BASE(3,:)=[2,1,3];
SEQ_ST.BASE(4,:)=[2,3,1];
SEQ_ST.BASE(5,:)=[3,1,2];
SEQ_ST.BASE(6,:)=[3,2,1];
SEQ_ST.NUMCHUNK=size(SEQ_ST.BASE,1);
%Ignoring chunk associations (12 sequences)
chunkpairs(1,:)=[1,2];chunkpairs(2,:)=[1,3];
chunkpairs(3,:)=[2,1];chunkpairs(4,:)=[2,5];
chunkpairs(5,:)=[3,4];chunkpairs(6,:)=[3,1];
chunkpairs(7,:)=[4,3];chunkpairs(8,:)=[4,6];
chunkpairs(9,:)=[5,6];chunkpairs(10,:)=[5,2];
chunkpairs(11,:)=[6,5];chunkpairs(12,:)=[6,4];
SEQ_ST.CHUNKPAIRS=chunkpairs;
SEQ_ST.NUMSEQTYPE=size(chunkpairs,1);
SEQ_ST.SERTYPE=[1:SEQ_ST.NUMSEQTYPE];%sequence type
SEQ_ST.MAXSEQLENGTH=[6];%smaximum elements in sequence
SEQ_ST.ONTRACKRATE=0.50;

baseseq=cell(size(SEQ_ST.CHUNKPAIRS,1),1);
for c=1:size(SEQ_ST.CHUNKPAIRS,1)
    c1=SEQ_ST.BASE(SEQ_ST.CHUNKPAIRS(c,1),:);
    c2=SEQ_ST.BASE(SEQ_ST.CHUNKPAIRS(c,2),:);
    baseseq(c)={horzcat({c1},{c2})};
end

SEQ_ST.BASESEQUENCE=vertcat(baseseq{:});

%%
% 1 experiment = 2 phases
% 1 phase = 4 runs 
% 1 run = 3 sequences (with 4 repetitions + 1)
% Make sure each element appears at each within-chunk position
clear seq
a=[SEQ_ST.BASESEQUENCE{4,:}]
b=[SEQ_ST.BASESEQUENCE{6,:}]
c=[SEQ_ST.BASESEQUENCE{12,:}]

seq=vertcat(a,b,c)
sum(seq,1)



