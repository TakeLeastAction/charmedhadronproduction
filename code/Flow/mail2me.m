function mail2me(subject,content)
MailAddress = 'offmsg@126.com';
password = 'QAZXSWEDC';
setpref('Internet','E_mail',MailAddress);
setpref('Internet','SMTP_Server','smtp.126.com');
setpref('Internet','SMTP_Username',MailAddress);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
sendmail('sunkaijia@sjtu.edu.cn',subject,content);
