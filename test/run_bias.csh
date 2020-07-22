#!/bin/csh -f

set ntoys = 10 # 5                                                                                                                                      
set ndirs = 500  # 2                                                                                                                                 
set index_fit = 0


foreach index_toy ( 1 2 )                                                                                                                                          
 #   echo "++++++ Produce toys with pdfindex = " ${index_toy} "++++++++"
    foreach mp ( 130 140 160 ) 

	foreach mu ( 0 1 2 3 )
	    mkdir SUB_BIAS_mu${mu}_M-${mp}_pdftoy_${index_toy}_pdffit_${index_fit}
	    cd SUB_BIAS_mu${mu}_M-${mp}_pdftoy_${index_toy}_pdffit_${index_fit}
	    echo Processing Mass point $mp GeV with injected signal strength  $mu
	    ../htc_all.csh $ntoys $ndirs $mu $mp $index_toy $index_fit
	    cd ..
	    echo "============================================="
	end
    end
end
