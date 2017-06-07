%import three variables from excel: Sequence1, CanonicalSequence1, 'Sequence1
%number'

%Sequence1=sequence;
%CanonicalSequence1=Precursor1;

%%
%%initialize output variables for the swalign() fucntion
Alignments = cell(length(CanonicalSequence1), 1);   
Starts = [];
index_empty=[];  %index of empty entries in 'Sequence1'---imported data
index_non_empty = [];


for i=1:length(CanonicalSequence1)
    if isempty(CanonicalSequence1{i})
        index_empty=[index_empty;i];
        Alignments{i} = '';
        Starts = [Starts; 0 0];
        continue;
    end
    display(i);
    
    
    %%perform local aligment
    index_non_empty = [index_non_empty;i];
    [~, Alignment, Start] = swalign(Sequence1{i}, CanonicalSequence1{i}, 'Alphabet', 'NT');
    Alignments{i} = Alignment;
    Starts = [Starts;Start'];      
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  process 5'  end%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First: 5'-overhang in Sequence1
my_starts = Starts(index_non_empty,:);
index_5 = my_starts(:, 1) > ones(length(index_non_empty),1) & my_starts(:, 2) == ones(length(index_non_empty),1);
index_5 = index_non_empty(index_5); %index in the original Sequence1 file
out_five_overhang = cell(length(Sequence1), 1);
for i=1:size(index_5, 1)
    
    start_position = Starts(index_5(i),1);  %temporaty
    seq = Sequence1{index_5(i)};                       %temporary
    %seq(1:start_position-1)
    out_five_overhang{index_5(i),1}=seq(1:start_position-1);
  
end

%% Second: 5'-end deletion
index_5_del = my_starts(:, 1) == ones(length(index_non_empty),1) & my_starts(:, 2) > ones(length(index_non_empty),1);
index_5_del = index_non_empty(index_5_del); %index in the original Sequence1 file
out_five_del = cell(length(Sequence1), 1);
for i=1:size(index_5_del, 1)
    
    start_position = Starts(index_5_del(i),2);  %temporaty
    seq = CanonicalSequence1{index_5_del(i)};                       %temporary
    %seq(1:start_position-1)
    out_five_del{index_5_del(i),1}=seq(1:start_position-1);
  
end

%% third: 5'-mutation
index_5_mut = my_starts(:, 1) > ones(length(index_non_empty),1) & my_starts(:, 2) > ones(length(index_non_empty),1);
index_5_mut = index_non_empty(index_5_mut); %index in the original Sequence1 file
out_five_mut_seq = cell(length(Sequence1), 1);
for i=1:size(index_5_mut, 1)
    
    start_position = Starts(index_5_mut(i),1);  %temporaty
    seq = Sequence1{index_5_mut(i)};                       %temporary
    %seq(1:start_position-1)
    out_five_mut_seq{index_5_mut(i),1}=seq(1:start_position-1);
  
end
out_five_mut_can = cell(length(Sequence1), 1);
for i=1:size(index_5_mut, 1)
    
    start_position = Starts(index_5_mut(i),2);  %temporaty
    seq = CanonicalSequence1{index_5_mut(i)};                       %temporary
    %seq(1:start_position-1)
    out_five_mut_can{index_5_mut(i),1}=seq(1:start_position-1);
  
end





%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% process middle %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
non_empty_Aligments = Alignments(index_non_empty,:);
%reshape the non-empty alignments
my_Alignments = cell(size(non_empty_Aligments,1), 3);

x=non_empty_Aligments(1);
y=non_empty_Aligments{1};

for i=1:size(non_empty_Aligments,1)
    current = non_empty_Aligments{i};
    my_Alignments{i,1}= current(1,:);
    my_Alignments{i,2}= current(2,:);
    my_Alignments{i,3}= current(3,:);
    
end
%% prepare alignment output (three columns)
out_alignments = cell(length(Sequence1),3);
for i=1:length(index_non_empty)
    out_alignments{index_non_empty(i),1}=my_Alignments{i,1};
    out_alignments{index_non_empty(i),2}=my_Alignments{i,2};
    out_alignments{index_non_empty(i),3}=my_Alignments{i,3};
    
end
%extract only the match symbols
 
match_symbols = my_Alignments(:,2);
idx_space = cellfun(@(x) sum(isspace(x)), match_symbols);
index_space=index_non_empty(logical(idx_space));
out_with_space = cell(length(Sequence1), 1);
for i=1:length(index_space)
    out_with_space{index_space(i),1}='yes';
end

%
idx_star = cellfun(@(x) sum(x=='.'), match_symbols);
idx_star = index_non_empty(logical(idx_star));
out_with_star = cell(length(Sequence1), 1);
for i=1:length(idx_star)
    out_with_star{idx_star(i),1}='yes';
end






%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% process 3' end %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First: 3'-overhang in Sequence1
%my_starts = Starts(index_non_empty,:);  %alread computed
%my_Alignments = cell(size(non_empty_Aligments,1), 3);   %already computed
%get length of alignments between 'Sequence1' and "CanonicalSequence1'
all_alignment_Sequence11=my_Alignments(:,1);
all_alignment_Sequence12=my_Alignments(:,3);
aligment_lengths1 = [];
aligment_lengths2 = [];

for i = 1:size(all_alignment_Sequence11, 1)
    temp_seq = all_alignment_Sequence11{i};
      
    aligment_lengths1=[aligment_lengths1;length(temp_seq)];
end


for i = 1:size(all_alignment_Sequence12, 1)
    temp_seq = all_alignment_Sequence12{i};
      
    aligment_lengths2=[aligment_lengths2;length(temp_seq)];
end
%sum(aligment_lengths1~=aligment_lengths2)
aligment_lengths=aligment_lengths1;


%% compute length of 'Sequence1' and 'CanonicalSequence1'
non_empty_Sequence1 = Sequence1(index_non_empty);
non_empty_can = CanonicalSequence1(index_non_empty);
length_non_empty_Sequence1=[];
length_non_empty_can=[];
for i=1:size(non_empty_Sequence1,1)
    length_non_empty_Sequence1=[length_non_empty_Sequence1;length(non_empty_Sequence1{i})];
end
for i=1:size(non_empty_can,1)
    length_non_empty_can=[length_non_empty_can;length(non_empty_can{i})];
end

%% 3'-overhang on 'sequnce'
index_Sequence1_3_over = my_starts(:, 1) + aligment_lengths -1 < length_non_empty_Sequence1 & my_starts(:, 2) + aligment_lengths -1 == length_non_empty_can;
index_Sequence1_3_over = index_non_empty(index_Sequence1_3_over); %index in the original Sequence1 file
out_Sequence1_three_overhang = cell(length(Sequence1), 1);
for i=1:size(index_Sequence1_3_over, 1)
    
    start_position = Starts(index_Sequence1_3_over(i),1) + size(Alignments{index_Sequence1_3_over(i),:},2);  %temporaty
    seq = Sequence1{index_Sequence1_3_over(i)};                       %temporary
    %seq(1:start_position-1)
    out_Sequence1_three_overhang{index_Sequence1_3_over(i),1}=seq(start_position:end);
    
  
end


%% 3'-overhang on canonical Sequence1
index_can_3_over = my_starts(:, 2) + aligment_lengths -1 < length_non_empty_can  & my_starts(:, 1) + aligment_lengths -1 == length_non_empty_Sequence1;
index_can_3_over = index_non_empty(index_can_3_over); %index in the original Sequence1 file
out_can_three_overhang = cell(length(Sequence1), 1);

for i=1:size(index_can_3_over, 1)
    
    start_position = Starts(index_can_3_over(i),1) + size(Alignments{index_can_3_over(i),:},2);  %temporaty
    seq = CanonicalSequence1{index_can_3_over(i)};                       %temporary
    %seq(1:start_position-1)
    out_can_three_overhang{index_can_3_over(i),1}=seq(start_position:end);
    
  
end

%% third: 3'-mutation
index_3_mut = my_starts(:, 1) + aligment_lengths -1 < length_non_empty_Sequence1 & my_starts(:, 2) + aligment_lengths -1 < length_non_empty_can;
index_3_mut = index_non_empty(index_3_mut); %index in the original Sequence1 file
out_three_mut_seq = cell(length(Sequence1), 1);
for i=1:size(index_3_mut, 1)
    
    
    
    start_position = Starts(index_3_mut(i),1) + size(Alignments{index_3_mut(i),:},2);  %temporaty
    seq = Sequence1{index_3_mut(i)};                       %temporary
    %seq(1:start_position-1)
    out_three_mut_seq{index_3_mut(i),1}=seq(start_position:end);
  
end
out_three_mut_can = cell(length(Sequence1), 1);
for i=1:size(index_3_mut, 1)
    
    
    
    start_position = Starts(index_3_mut(i),2) + size(Alignments{index_3_mut(i),:},2);  %temporaty
    seq = CanonicalSequence1{index_3_mut(i)};                       %temporary
    %seq(1:start_position-1)
    out_three_mut_can{index_3_mut(i),1}=seq(start_position:end);
  
end


%%
%combine all the output
 final_output = [out_alignments, out_five_overhang, out_five_del, out_five_mut_seq, out_five_mut_can, out_Sequence1_three_overhang,out_can_three_overhang, out_three_mut_seq, out_three_mut_can, out_with_space];

 
 
%% save the final_output as a textfile
T = cell2table(final_output,'VariableNames',{'out_alignments1', 'out_alignments2', 'out_alignments3','out_five_overhang', 'out_five_del', 'out_five_mut_seq', 'out_five_mut_can', 'out_Sequence1_three_overhang','out_can_three_overhang', 'out_three_mut_seq', 'out_three_mut_can','out_with_space'});

writetable(T,'B2Precursor1_new.dat')

%%









