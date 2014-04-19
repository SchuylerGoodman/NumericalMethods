function num = trunc(x, y)
%trunc
%truncates a number x to a a maximum of y digits
time1=clock;
num=x;
dig=0;
first=0;
neg=x<0;

while fix(num)~=0
    num=fix(num/10);
    if first==1
        dig=dig+1;
    end
    first=1;
end
if dig==0
    num=x;
    while fix(num)==0
        num=num*10;
        dig=dig-1;
    end
end
digits=['%.' num2str(y+abs(dig)) 'f'];
numst=num2str(x,digits);
dec=dig<0;
place=1+neg;
if dec==1
    place=2+abs(dig)+neg;
end
pre=numst(place);
if neg==1
    pre=strcat(numst(1),pre);
end
post='.';
count=1;
if length(numst)>(place)
    for i=(place+1):1:length(numst)
        if count>=y
            break;
        end
        if numst(i)~='.'
            post=strcat(post,numst(i));
            count=count+1;
        end
    end
end
num=str2double(strcat(pre,post,'e',num2str(dig)));
time2=clock;
%fprintf('elapsed time: %g seconds\n', etime(time2,time1));
end

