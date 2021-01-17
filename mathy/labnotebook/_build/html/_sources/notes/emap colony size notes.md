# 2020-03-27 EMAP Colony Size Notes

- ReplicateControlCorrWithRowMut() function gives correlation coefficients.
    - Use the code from `ReplicateControlCorrWithRowMut.m`v to find how to get raw colony sizes to compute yourself.
- In the code Tina ran to compute the replicate control correlation she ran the following in MATLAB:

``` MATLAB
load('/Users/tina/Documents/GSP1/E_MAP_data/12-25-2014/December2014_screen.mat');
ReplicateControlCorrWithRowMut(rawN);
```

- So she loaded in the matrix, and then passed rawN through the function:

``` MATLAB
function ReplicateControlCorrWithRowMut(rawN_int)
%This fcn computes the correlation coefficients between the colony sizes on
%all plates with the median colony size (controlsize) in each corresponding 
%position. The important thing here is that all 3 replicates should give 
%similar correlation coefficients. Also, if they are lower than say 0.2,
%something may be off. If corr coeff is e.g. 0.3 but all three similar, 
%this is probably just a strain with a lot of interactions.

%originally by Sean, slightly adapted by Amelie (2013)


for i=1:1:length(rawN_int.rowlabels)
    fprintf ('%i\t%s-%s',i,char(getGeneName(rawN_int,i)), rawN_int.rowMut{i});
    for j=1:1:size(rawN_int.size,3)
        fprintf ('\t%.3f', myNanSpearman(rawN_int.size(i,:,j),rawN_int.controlsize));
    end
    fprintf ('\n');
end
```

``` python
def read_in_EMAP_screen(filename)

    mdict = {}
    mat = scipy.io.loadmat(file_name=filename, mdict=mdict)

    # query ORF, 33 x 1
    rowlabels = [label[0][0] for label in mdict['rawN']['rowlabels'][0][0]]

    # query mutation type (including Gsp1 point mutant labels), 33 x 1
    rowMut = [label[0][0] for label in mdict['rawN']['rowMut'][0][0]]

    # library ORF, 1536 x 1
    collabels = [label[0][0] for label in mdict['rawN']['collabels'][0][0]]

    # library mutation type, 1536 x 1
    colMut = [label[0][0] for label in mdict['rawN']['colMut'][0][0]]

#     # colony sizes, 33 x 1536 x 3
#     colony_size_matrix = mdict['rawN']['size'][0][0]

#     rep1 = pd.DataFrame(colony_size_matrix[:,:,0])
#     rep2 = pd.DataFrame(colony_size_matrix[:,:,1])
#     rep3 = pd.DataFrame(colony_size_matrix[:,:,2])


    # average colony size, 33 x 1536
    avg_colony_size = mdict['rawN']['meansize'][0,0]

    # control colony size, 33 x 1536
    ctrl_colony_size = mdict['rawN']['controlsize'][0,0]

    # average colony size, 33 x 1536
    sd_colony_size = mdict['rawN']['sdsize'][0,0]

    # average colony size, 33 x 1536
    sd_ctrl_colony_size = mdict['rawN']['controlsd'][0,0]


    ORF2gene = dict(zip([ORF[0][0] for ORF in mdict['rawN']['geneToOrf'][0,0]['orfname'][0,0]],
                        [ORF[0][0] for ORF in mdict['rawN']['geneToOrf'][0,0]['genename'][0,0]]))
```