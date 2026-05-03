load factors

%load VOI_V2_1
load VOI_V1_1
v1=Y;

load VOI_V5_1
v5=Y;

subplot(3,1,1);
plot(v1);
set(gca,'XTickLabel',[]);
%title('V2');
title('V1');

subplot(3,1,2);
plot(v5);
set(gca,'XTickLabel',[]);
title('V5');

subplot(3,1,3);
plot(attention);


title('Attention');
axis([0 400 -0.5 1.5]);

no_attention=abs(attention-motion);


v1_att=attention.*v1;
v5_att=attention.*v5;
v1_att=v1_att(find(v1_att~=0));
v5_att=v5_att(find(v5_att~=0));

v1_natt=no_attention.*v1;
v5_natt=no_attention.*v5;
v1_natt=v1_natt(find(v1_natt~=0));
v5_natt=v5_natt(find(v5_natt~=0));

cc=corrcoef([v1_natt,v5_natt]);
c=cc(1,2);
disp(sprintf('Correlation when not attending is %1.3f',c));

cc=corrcoef([v1_att,v5_att]);
c=cc(1,2);
disp(sprintf('Correlation when attending is %1.3f',c));


