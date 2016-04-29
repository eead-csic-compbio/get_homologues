# LaTeX2HTML 2008 (1.71)
# Associate images original text with physical files.


$key = q/G;MSF=1.6/;
$cached_env_img{$key} = q|<IMG
 WIDTH="16" HEIGHT="29" ALIGN="MIDDLE" BORDER="0"
 SRC="|."$dir".q|img1.png"
 ALT="$ G$">|; 

$key = q/G^{2};MSF=1.6/;
$cached_env_img{$key} = q|<IMG
 WIDTH="22" HEIGHT="35" ALIGN="MIDDLE" BORDER="0"
 SRC="|."$dir".q|img2.png"
 ALT="$ G^{2}$">|; 

$key = q/Evalue<max(Evalue);MSF=1.6/;
$cached_env_img{$key} = q|<IMG
 WIDTH="155" HEIGHT="31" ALIGN="MIDDLE" BORDER="0"
 SRC="|."$dir".q|img3.png"
 ALT="$ Evalue &lt; max(Evalue)$">|; 

1;

