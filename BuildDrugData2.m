%Eine Tropennacht ist eine Nacht (18 bis 06 UTC), in der das Minimum der Lufttemperatur >= 20 C beträgt.
%https://www.dwd.de/DE/klimaumwelt/klimaatlas/klimaatlas_node.html
%eventuell auch hier heise tage prognose und daten her

%https://www.umweltbundesamt.de/indikator-heisse-tage#Die wichtigsten Fakten
%hitzetag Anzahl der Tage mit einem Lufttemperatur-Maximum über 30 Grad Celsius
%(Gebietsmittel) in deutschland 1951-2017


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
%interpolate zig
%1964-2017

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