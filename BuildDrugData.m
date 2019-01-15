%Eine Tropennacht ist eine Nacht (18 bis 06 UTC), in der das Minimum der Lufttemperatur >= 20 C beträgt.
%https://www.dwd.de/DE/klimaumwelt/klimaatlas/klimaatlas_node.html
%eventuell auch hier heise tage prognose und daten her

%https://www.umweltbundesamt.de/indikator-heisse-tage#Die wichtigsten Fakten
%hitzetag Anzahl der Tage mit einem Lufttemperatur-Maximum über 30 Grad Celsius
%(Gebietsmittel) in deutschland 1951-2017
heatdays=[3.02
7.91
5.08
2.53
0.93
0.58
6.85
1.58
5.63
1.3
3.97
2.34
4.5
9.93
1.18
2.64
3.95
2.38
6.63
2.01
7.08
4.7
5.78
2.75
6.63
10.2
1.25
1.89
2.29
1.45
2.19
6.54
9.93
3.3
2.63
4.31
1.64
2.1
4.67
5.91
5.27
9.64
1.89
16.27
10.54
3
5.16
7.08
5.16
5.67
7.36
6.23
19.01
4.7
7.57
14.34
4.82
7.24
4.58
10.63
4.27
7.62
10.47
6.1
17.6
9.2
6.8
];

%1951-2018 data from email dwd
heatnights=[0.01
0.75
0.01
0.01
0
0
0.63
0
0.32
0
0
0.01
0.15
0.27
0.01
0.04
0.07
0.03
0.43
0.01
0.2
0.2
0.1
0.08
0.22
0.51
0
0.03
0.02
0.02
0.03
0.22
0.66
0.29
0.04
0.22
0.19
0.12
0.17
0.12
0.48
0.64
0
1.7
0.19
0.01
0.21
0.25
0.07
0.1
0.28
0.28
1.54
0.08
0.17
1
0.18
0.03
0.01
0.84
0.06
0.2
0.5
0.22
1.3
0.08
0.04
1.34];

%https://de.statista.com/statistik/daten/studie/4917/umfrage/inflationsrate-in-deutschland-seit-1948/
%inflationsrate 1950-2017 westdeutschland 
%Veränderung des Verbraucherpreisindex im Vergleich zum Vorjahr
inflation=[1.8
0.5
0.3
0.9
1.5
2
2.1
1.1
0.3
2.6
2.3
1.5
1.6
1.6
1.1
1.4
2
1.4
0.6
1
2
1.4
1.8
2.6
4.5
5.1
3.7
2.6
2.8
1.2
0.2
-0.1
2
2.5
3.2
5.2
6.3
5.4
4.1
2.7
3.7
4.2
6
6.9
7.1
5.4
5.2
3.6
1.8
1.6
1.9
3.3
3.2
2.4
3
2.8
2.5
1.6
0.6
2.3
2
2.8
1.4
0.4
-1.7
2.1
7.6
-6.4
]';
inflation=fliplr(inflation);
inflation=inflation';
%https://de.statista.com/statistik/daten/studie/5382/umfrage/alkoholverbrauch-je-einwohner-an-reinem-alkohol/
%average alkoholverbrauch je einwohner ab 15 jahren in liter im jahr
%ganz deutschland
%1970,1980,1990,2000,2010,2011,2013,2015,2017
alk=[14.4,15.1,13.4,12,10.7,11,10.7,10.7,10.7]';
hh=[1970,1980,1990,2000,2010,2011,2013,2015,2017];
%plot(hh,alk,'g');
hold on
%https://www.gesis.org/fileadmin/upload/dienstleistung/daten/soz_indikatoren/Schluesselindikatoren/G028.pdf
%1950-1989 west germany 1990-2005 germany
%1. Angenommener Alkoholgehalt fr Bier: 4%, fr Wein und Schaumwein: 10%, fr Branntwein: direkte Ausweisung der Alkoholmenge.                                          
%2. In die Berechnung eingegangene Werte fr Bier erst ab 1993 ohne alkoholfreies Bier.   
%3. West erst ab 1960 incl. Saarland; 
%4. West erst ab 1952 inkl. Berlin-West. 
%5 1990 missing, take average of 1989 and 1991
%daily pure alkoholkonsum in ml of average person
alk_2=[ 12.1 ,13.3, 13.0,13.8 ,15.0,15.9,16.9,19.2,19.9,21.2,23.6,26.2,27.8,29.4,29.9,32.0,31.6,31.4,33.0,34.6,36.4,38.7,37.7,38.8,36.8,39.3,40.5,38.5,38.2,38.9,38.2,37.0,35.9,36.0,34.6,34.6,33.6,33.9,33.3,32.9,-1,34.6,33.1,32.6,32.6,32.0,31.7,31.4,30.9,30.6,30.7,30.2,29.6,28.9,28.5,28.1 ];
alk_2(41)=(alk_2(40)+alk_2(39))/2;
%convert to yearly liter
%alk_2=alk_2*365/1000;


%https://de.statista.com/statistik/daten/studie/5384/umfrage/verbrauch-je-einwohner-an-alkohol-in-deutschland-seit-1990/
%weiterfuhrung 2008-2017 daten konsum in liter pro jahr pro kopf
%annahme spirituosen 30%
%bier=[111.1,109.9,107.4,109.3,107.3,106.6,106.9,105.9,104.1,101.1]*0.04;
%wein =[20.7,20.1,20.5,20.6,20.8,21.1,20.7,20.5,21.1,20.9]*0.1;
%schaumwein=[3.9,3.9,3.9,4.2,4.2,4,3.9,3.7,3.7,3.9]*0.1;
%schnaps=[5.5,5.4,5.4,5.5,5.5,5.5,5.4,5.4,5.4,5.4]*0.4;
%alk_est=bier+wein+schaumwein+schnaps;
%hhh=linspace(2008,2017,10);
%plot(hhh,alk_est);

%temp=linspace(alk_2(end),alk_est(1),4);
%alk_2=[alk_2, temp(2),temp(3),alk_est];
%hhhh=linspace(1950,2017,68);
%plot(hhhh,alk_2,':')
%alk_2=alk_2';
%https://de.statista.com/statistik/daten/studie/444502/umfrage/raucheranteil-unter-jugendlichen-und-jungen-erwachsenen-nach-geschlecht/
%Raucheranteil 12-25 jahre
%1979,1982,1986,1989,1993,1997,2001,2004,2008,2010,2011,2012,2014,2015
%ab 1993 mit neue bundeslaender 
rauch_m_young=[47.3,42.2,45.8,46.4,40.4,42.8,38.3,36,31.2,30.5,26.9,27.4,24.5,20.1]';
rauch_f_young=[40.2,38.9,44.3,39.9,33.8,39.4,36.1,35,32.8,26.2,27.3,24.9,19.4,17.7]';

%https://de.statista.com/statistik/daten/studie/444495/umfrage/raucheranteil-unter-jugendlichen-und-jungen-erwachsenen-in-deutschland/
%ab 1993 mit neue bundeslaender 
rauch_ges_young=[43.9,40.6,45.1,43.3,37.2,41.2,37.2,35.5,32,28.4,27.1,26.2,22,18.9]';

%raucheranteil ab 14 jahre westdeutsche bevolkerung
%https://www.ifd-allensbach.de/uploads/tx_reportsndocs/prd_0325.pdf
%1950 1955 1965 1975 1985 1996 2000 2001 2002 2003
rauch_m=[88,83,74,60,47,39,41,40,40,39]';
rauch_f=[21,21,24,29,28,27,27,27,28,28]';
rauch_ges=[51,49,46,44,37,32,33,33,34,33]';
x_org=[1950,1955,1965,1975,1985,1996,2000,2001,2002,2003];
rauch_m_org=[88,83,74,60,47,39,41,40,40,39]';
rauch_f_org=[21,21,24,29,28,27,27,27,28,28]';
rauch_ges_org=[51,49,46,44,37,32,33,33,34,33]';


%https://de.statista.com/statistik/daten/studie/6187/umfrage/absatz-von-versteuerten-zigaretten-seit-1964/#0
%absatz verstuererter zigaretten in deutschland in mio stueck
%1964 1970 1975 1980 1985 1991 1995 2000 2005 2006 2007 2008 2009 2010 2011
%2012 2013 2014 2015 2016 2017
zig=[90156,118051,126200,128425,120408,146480,135029,139625,95827,93465,91497,87979,86607,83565,87556,82405,80266,79521,81267,75016,75838]';

%i will use linear interpolation in beetwen years to fill out the data
%interpolate alk
%1970-2017
h1=linspace(14.4,15.1,11);
h2=linspace(15.1,13.4,11);
h2=h2(2:end);
h3=linspace(13.4,12,11);
h3=h3(2:end);
h4=linspace(12,10.7,11);
h4=h4(2:end);
alk=[h1,h2,h3,h4,11,(11+10.7)/2,10.7,10.7,10.7,10.7,10.7]';

%interpolate zig
%1964-2017
h1=linspace(90156,118051,7);
h2=linspace(118051,126200,6);
h2=h2(2:end);
h3=linspace(126200,128425,6);
h3=h3(2:end);
h4=linspace(128425,120408,6);
h4=h4(2:end);
h5=linspace(120408,146480,7);
h5=h5(2:end);
h6=linspace(146480,135029,5);
h6=h6(2:end);
h7=linspace(135029,139625,6);
h7=h7(2:end);
h8=linspace(139625,95827,6);
h8=h8(2:end);
h9=zig(10:end)';
zig=[h1,h2,h3,h4,h5,h6,h7,h8,h9]';

%interpolate rauch_young
%1979-2015
h1=linspace(47.3,42.2,4);
h2=linspace(42.2,45.8,5);
h2=h2(2:end);
h3=linspace(45.8,46.4,4);
h3=h3(2:end);
h4=linspace(46.4,40.4,5);
h4=h4(2:end);
h5=linspace(40.4,42.8,5);
h5=h5(2:end);
h6=linspace(42.8,38.3,5);
h6=h6(2:end);
h7=linspace(38.3,36,4);
h7=h7(2:end);
h8=linspace(36,31.2,5);
h8=h8(2:end);
h9=linspace(31.2,30.5,3);
h9=h9(2:end);
rauch_m_young=[h1,h2,h3,h4,h5,h6,h7,h8,h9,26.9,27.4,(27.4+24.5)/2,24.5,20.1]';

h1=linspace(40.2,38.9,4);
h2=linspace(38.9,44.3,5);
h2=h2(2:end);
h3=linspace(44.3,39.9,4);
h3=h3(2:end);
h4=linspace(39.9,33.8,5);
h4=h4(2:end);
h5=linspace(33.8,39.4,5);
h5=h5(2:end);
h6=linspace(39.4,36.1,5);
h6=h6(2:end);
h7=linspace(36.1,35,4);
h7=h7(2:end);
h8=linspace(35,32.8,5);
h8=h8(2:end);
rauch_f_young=[h1,h2,h3,h4,h5,h6,h7,h8,(32.8+26.2)/2,26.2,27.3,24.9,(24.9+19.4)/2,19.4,17.7]';

h1=linspace(43.9,40.6,4);
h2=linspace(40.6,45.1,5);
h2=h2(2:end);
h3=linspace(45.1,43.3,4);
h3=h3(2:end);
h4=linspace(43.3,37.2,5);
h4=h4(2:end);
h5=linspace(37.2,41.2,5);
h5=h5(2:end);
h6=linspace(41.2,37.2,5);
h6=h6(2:end);
h7=linspace(37.2,35.5,4);
h7=h7(2:end);
h8=linspace(35.5,32,5);
h8=h8(2:end);
rauch_ges_young=[h1,h2,h3,h4,h5,h6,h7,h8,(28.4+32)/2,28.4,27.1,26.2,(22+26.2)/2,22,18.9]';

%interpolate rauch 
%1950-2003
h1=linspace(88,83,6);
h2=linspace(83,74,11);
h2=h2(2:end);
h3=linspace(74,60,11);
h3=h3(2:end);
h4=linspace(60,47,11);
h4=h4(2:end);
h5=linspace(47,39,12);
h5=h5(2:end);
h6=linspace(39,41,5);
h6=h6(2:end);
rauch_m=[h1,h2,h3,h4,h5,h6,40,40,39]';

h1=linspace(21,21,6);
h2=linspace(21,24,11);
h2=h2(2:end);
h3=linspace(24,29,11);
h3=h3(2:end);
h4=linspace(29,28,11);
h4=h4(2:end);
h5=linspace(28,27,12);
h5=h5(2:end);
h6=linspace(27,27,5);
h6=h6(2:end);
rauch_f=[h1,h2,h3,h4,h5,h6,27,28,28]';

h1=linspace(51,49,6);
h2=linspace(49,46,11);
h2=h2(2:end);
h3=linspace(46,44,11);
h3=h3(2:end);
h4=linspace(44,37,11);
h4=h4(2:end);
h5=linspace(37,32,12);
h5=h5(2:end);
h6=linspace(32,33,5);
h6=h6(2:end);
rauch_ges=[h1,h2,h3,h4,h5,h6,33,34,33]';

clearvars h1 h2 h3 h4 h5 h6 h7 h8 h9 bier wein schnaps schaumwein h hh hhh hhhh alk_est temp