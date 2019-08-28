function report2me(text)

url='http://2.smsfx.sinaapp.com/send.php?';
tel='tel=13472616605';
pwd='pwd=667667xn';
aim='aim=13472616605';
%text='Be yourself, everyone else is already taken.?';
eml=strcat(url,tel,'&',pwd,'&',aim,'&','text=',text);
web(eml)