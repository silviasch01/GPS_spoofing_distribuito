1;
function thr = upper_quantile(x,Pfa), x=sort(x(:)); n=numel(x); thr=x(min(n,max(1,ceil((1-Pfa)*n)))); end
function thr = lower_quantile(x,Pfa), x=sort(x(:)); n=numel(x); thr=x(min(n,max(1,floor(Pfa*n)+1))); end
cells={'g10','alt500','alt1000','alt2000','alt3500','alt5000','alt7000','alt10000','loit_circle','loit_los2000','loit_los5000','loit_los10000'};
fid=fopen('permission_pd.csv','w'); fprintf(fid,'cell,mission,pd5,pd3,stat5_med,ratio_med\n');
for c=cells
  L=load(sprintf('cell_%s.mat',c{1}));
  if exist(sprintf('cell_%s_b.mat',c{1}),'file'), B=load(sprintf('cell_%s_b.mat',c{1})); for q=1:2, L.T{1,q}=[L.T{1,q};B.T{1,q}]; L.RID{1,q}=[L.RID{1,q};B.RID{1,q}]; end, end
  T0=L.T{1,1}; T1=L.T{1,2}; if size(T0,2)<5, T0(:,5)=abs(T0(:,3)-1); T1(:,5)=abs(T1(:,3)-1); end
  thr5=upper_quantile(T0(:,5),0.01); thr3=lower_quantile(T0(:,3),0.01);
  for r=1:20
    ii=L.RID{1,2}==r; if ~any(ii), continue; end
    fprintf(fid,'%s,%d,%.4f,%.4f,%.4f,%.4f\n',c{1},r,mean(T1(ii,5)>thr5),mean(T1(ii,3)<thr3),median(T1(ii,5)),median(T1(ii,3)));
  end
end
fclose(fid); disp('PERMISSION DONE');
