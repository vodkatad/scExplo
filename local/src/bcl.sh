#bcl2fastq --sample-sheet /mnt/cold1/upload/scRNA_Candiolo_12_2021/ll/211209_A00959_0109_BHL7L3DRXY/SampleSheet.csv -o test -R /mnt/cold1/upload/scRNA_Candiolo_12_2021/ll/211209_A00959_0109_BHL7L3DRXY/
bcl2fastq --use-bases-mask=Y28,I10,I10,Y90 \
	  --create-fastq-for-index-reads \
	    --minimum-trimmed-read-length=8 \
	      --mask-short-adapter-reads=8 \
	        --ignore-missing-positions \
		  --ignore-missing-controls \
		    --ignore-missing-filter \
		      --ignore-missing-bcls \
		        -r 6 -w 6 \
			  -R /mnt/cold1/upload/scRNA_Candiolo_12_2021/ll/211209_A00959_0109_BHL7L3DRXY/ \
			    --output-dir=test \
			      --interop-dir=/mnt/cold1/upload/scRNA_Candiolo_12_2021/ll/211209_A00959_0109_BHL7L3DRXY/InterOp \
			        --sample-sheet=/mnt/cold1/upload/scRNA_Candiolo_12_2021/ll/211209_A00959_0109_BHL7L3DRXY/SampleSheet.csv
  
