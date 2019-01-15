%https://www.destatis.de/DE/ZahlenFakten/GesamtwirtschaftUmwelt/VGR/Inlandsprodukt/Tabellen/Volkseinkommen1925_pdf.pdf?__blob=publicationFile
%non optimal solution, very problematic, this data needs help
%actuallt GDP is not used, instead the related Volkseinkommen

%g_full= Volkseinkommen in jeweiligen Preisen je Einwohner (all years)

%1950-1960, westdeutschland ohne saarland und westberlin
g_full=[856,1028,1160,1230,1301,1477,1621,1754,1855,1996,2226];

%1960-1970 west deutschland
%1960 different value 2215 than above. 
%replace to match mortality data which is west geermany based
g_full(end)=2215;

g_full=[g_full,2373,2538,2653,2886,3127,3283,3284,3593,3941,4471];

%1970-1991 west deutschland values of VGR-Revision 2005
%1970 value of 4652 greatly different. I will use the average of both
g_full(end)=(4471+4652)/2;

g_full=[g_full,5080,5507,6136,6630,6962,7602,8097,8667,9343,9897,10302,10701,11226,11876,12493,13194,13514,14290,15125,16092,16996];

%1992-2017 only data for the hole of germany is available 
%better than nothing for now... 1991 value is 15337
%g_full=log([g_full,16231,16391,16951,17577,17787,18057,18432,18671,19089,19589,19695,19776,20779,21108,22312,23241,23487,22633,23955,25264,25546,26089,26935,27727,28391,29448]');

%found the data so above not the case anymore
%einwohner westdeutschland 1991-2017 in 1000
einwohner=[ 65348.7  
 66066.3  
 66628.7  
 66921.5  
 67155.9  
 67376.0  
 67477.3  
 67483.5  
 67539.5  
 67668.6  
 67850.8  
 68040.5  
 68125.1  
 68136.2  
 68117.4  
 68059.3  
 67993.0  
 67882.8  
 67712.9  
 67607.4  
 67671.2  
 67874.6  
 68131.9  
 68481.5  
 69135.2  
 69758.8  
 70096.2  
];

%volkseinkommen westdeutchland 1991-2016 in mio euro

Neinkommen=[   17216  
  17941  
  17878  
  18299  
  18889  
  19066  
  19352  
  19741  
  19963  
  20406  
  20892  
  20961  
  21008  
  22036  
  22412  
  23666  
  24598  
  24799  
  23808  
  25202  
  26557  
  26830  
  27362  
  28243  
  29046  
  29719  
];

%17216 vs 16996 for 1991. i choose average
g_full(end)=(g_full(end)+Neinkommen(1))/2;
g_full=([g_full,Neinkommen(2:end)']');

%cut out the years used in the Lee Carter mortality data and log
g_t=log (g_full(t_min-1950+1:t_min-1950+t));


%g_t=(1+r_t)g_t-1 from 1951 to 2016
r_full=zeros(length(g_full)-1,1);
for i=1:length(r_full)
    r_full(i)= g_full(i+1)/g_full(i)-1;
end
g_full=log(g_full);   
    

clearvars Neinkommen einwohner g_calc
