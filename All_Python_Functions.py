def mergePeaks(peakFolder, macs2FilePattern):
    mergedPeakFile = peakFolder+'mergedPeaks.bed'
    tmpPeakFile = peakFolder + 'tmpPeak.txt'
    concatCmd='cat '
    tmpPeakSortFile = peakFolder + 'tmpPeakSort.txt'
    for root, folder, files in os.walk(peakFolder):
        files = [os.path.join(root, f) for f in files if (f.endswith(macs2FilePattern))]
        for f in files:
            concatCmd = concatCmd + f + ' '

    concatCmd = concatCmd + ' > ' + mergedPeakFile
    subprocess.call(concatCmd, shell=True)
    sortCmd = 'sort -k1,1 -k2,2n ' + mergedPeakFile + ' > ' + tmpPeakFile  
    subprocess.call(sortCmd, shell=True)
    if(macs2FilePattern == '_peaks.narrowPeak' or macs2FilePattern == 'peaks.narrowPeak'):
        bedToolsCmd ='bedtools merge -i ' + tmpPeakFile + ' -c 9 -o median > ' + mergedPeakFile
        subprocess.call(bedToolsCmd, shell=True)
    elif(macs2FilePattern == '_summits_shifted.bed'):
        bedToolsCmd ='bedtools merge -i ' + tmpPeakFile + ' -c 5 -o median > ' + mergedPeakFile
        subprocess.call(bedToolsCmd, shell=True)
    else:
        print('Invalid pattern... Allowed only _summits_shifted.bed or _peaks.narrowPeak... Exiting!!')
        return
    sortMergeCmd = 'sort -nrk4 ' +  mergedPeakFile + ' > ' + tmpPeakSortFile
    subprocess.call(sortMergeCmd, shell=True)
    #copyTopCmd = 'head -n ' + str(topX) + ' ' + tmpPeakSortFile + ' > ' + mergedPeakFile
    copyTopCmd =  'cat ' + tmpPeakSortFile + ' > ' + mergedPeakFile
    subprocess.call(copyTopCmd, shell=True)
    rmvCmd = 'rm ' + tmpPeakFile
    subprocess.call(rmvCmd, shell=True)
    rmvCmd = 'rm ' + tmpPeakSortFile
    subprocess.call(rmvCmd, shell=True)  

def getMergedBam(BAMFolder, BAMFilePattern):
    
    import pandas as pd
    
    mergeFiles = ''
    for root, folder, files in os.walk(BAMFolder):
            files = [os.path.join(root, f) for f in files if (f.endswith(BAMFilePattern))]
            for f in files:
                mergeFiles = mergeFiles +  f + ' '
    
        
    mergeCommand = 'samtools merge ' + BAMFolder + 'AllFiles_Filtered_Merged_nodup.bam' + ' ' + mergeFiles    
    subprocess.call(mergeCommand, shell=True)   
    
    sortMergeCmd = 'samtools sort ' + BAMFolder + 'AllFiles_Filtered_Merged_nodup.bam' \
                    + ' ' + BAMFolder + 'AllFiles_Filtered_Merged_nodup_sorted'        
    subprocess.call(sortMergeCmd, shell=True)
    
    indexMergeCmd = 'samtools index ' + BAMFolder + 'AllFiles_Filtered_Merged_nodup_sorted.bam'
    subprocess.call(indexMergeCmd, shell=True)


def getAggregatedPeak(BAMFolder, BAMFilePattern):
    getMergedBam(BAMFolder, BAMFilePattern)
    
    MergedMacs2Cmd = 'macs2 callpeak -t ' + BAMFolder + 'AllFiles_Filtered_Merged_nodup_sorted.bam' +' -n ' \
                    + BAMFolder + 'AllFiles_Filtered' +' -q 0.01 -g hs -f BAM --nomodel --nolambda --shift -75 \
                    --extsize 150 -B --SPMR --call-summits --keep-dup all'
    subprocess.call(MergedMacs2Cmd, shell=True)
