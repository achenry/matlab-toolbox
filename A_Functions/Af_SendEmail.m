function Af_SendEmail(recipient,subject,message,varargin)
% Function will send email from daniel.zalkind@colorado.edu (for now) to the (string)
% recipient with the (string) message and (strings to) file attachments

%% Set Email Preferences
% Server, Username and Password (make this no read besides me)

mail = 'dan.zalkind@gmail.com';
psswd = '1bladedWTs!';
host = 'smtp.gmail.com';
port  = '465';

setpref( 'Internet','E_mail', mail );
setpref( 'Internet', 'SMTP_Server', host );
setpref( 'Internet', 'SMTP_Username', mail );
setpref( 'Internet', 'SMTP_Password', psswd );

%Java properties
props = java.lang.System.getProperties;
props.setProperty( 'mail.smtp.user', mail );
props.setProperty( 'mail.smtp.host', host );
props.setProperty( 'mail.smtp.port', port );
props.setProperty( 'mail.smtp.starttls.enable', 'true' );
props.setProperty( 'mail.smtp.debug', 'true' );
props.setProperty( 'mail.smtp.auth', 'true' );
props.setProperty( 'mail.smtp.socketFactory.port', port );
props.setProperty( 'mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory' );
props.setProperty( 'mail.smtp.socketFactory.fallback', 'false' );

%% Send Email (no attachments yet)
if isempty(varargin)
    sendmail(recipient,subject,message);
else
    sendmail(recipient,subject,message,varargin);
end

