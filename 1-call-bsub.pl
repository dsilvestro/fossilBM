  #!/usr/bin/env perl

use File::Find;


        $indexVar2=1;
        while ($indexVar2<=1)
        {

$indexVar1=3;
while ($indexVar1<=3 )
{

	$variable="#!/bin/bash \n#BSUB -L /bin/bash\n#BSUB -cwd /scratch/ul/monthly/nsalamin/fossBM/ -n 1 -e err%J.txt -o out%J.txt\nmodule add R/3.0.2 \nRscript test.R  \n"; 
	system("echo \"$variable\" | bsub -q normal -n 1");



        $indexVar1++;
        }
 

	$indexVar2++;
        
}	
exit(0);

        
# "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o /Volumes/DSILVESTRO/data/mcmc_BM/0shift.pdf *_x_1.pdf
# "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o /Volumes/DSILVESTRO/data/mcmc_BM/1shift.pdf *_x_10.pdf
# "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o /Volumes/DSILVESTRO/data/mcmc_BM/2shift.pdf *_x8.pdf