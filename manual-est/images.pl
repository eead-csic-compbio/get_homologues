# LaTeX2HTML 2008 (1.71)
# Associate images original text with physical files.


$key = q/<95%;MSF=1.6/;
$cached_env_img{$key} = q|<IMG
 WIDTH="50" HEIGHT="30" ALIGN="MIDDLE" BORDER="0"
 SRC="|."$dir".q|img2.png"
 ALT="$ &lt;95\%$">|; 

$key = q/>75%;MSF=1.6/;
$cached_env_img{$key} = q|<IMG
 WIDTH="50" HEIGHT="30" ALIGN="MIDDLE" BORDER="0"
 SRC="|."$dir".q|img3.png"
 ALT="$ &gt;75\%$">|; 

$key = q/G;MSF=1.6/;
$cached_env_img{$key} = q|<IMG
 WIDTH="16" HEIGHT="31" ALIGN="MIDDLE" BORDER="0"
 SRC="|."$dir".q|img4.png"
 ALT="$ G$">|; 

$key = q/G^{2};MSF=1.6/;
$cached_env_img{$key} = q|<IMG
 WIDTH="23" HEIGHT="37" ALIGN="MIDDLE" BORDER="0"
 SRC="|."$dir".q|img5.png"
 ALT="$ G^{2}$">|; 

$key = q/geq40;MSF=1.6/;
$cached_env_img{$key} = q|<IMG
 WIDTH="36" HEIGHT="31" ALIGN="MIDDLE" BORDER="0"
 SRC="|."$dir".q|img1.png"
 ALT="$ \geq 40$">|; 

$key = q/.slashinstall.pl;MSF=1.6/;
$cached_env_img{$key} = q|<IMG
 WIDTH="78" HEIGHT="32" ALIGN="MIDDLE" BORDER="0"
 SRC="|."$dir".q|img6.png"
 ALT="$ ./install.pl$">|; 

1;

