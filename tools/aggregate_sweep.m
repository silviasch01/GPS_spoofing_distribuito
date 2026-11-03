1;
% aggregate the altitude-sweep and loitering/LOS cells: Pd (+CI), AUC for
% the 6 detectors, median ratio and distortion; writes sweep5.csv
function thr = upper_quantile(x,Pfa), x=sort(x(:)); n=numel(x); thr=x(min(n,max(1,ceil((1-Pfa)*n)))); end
function thr = lower_quantile(x,Pfa), x=sort(x(:)); n=numel(x); thr=x(min(n,max(1,floor(Pfa*n)+1))); end
function pd = pd_at_pfa(t0,t1,Pfa,small)
  t0=t0(~isnan(t0)); t1=t1(~isnan(t1));
  if small, thr=lower_quantile(t0,Pfa); pd=mean(t1<thr); else, thr=upper_quantile(t0,Pfa); pd=mean(t1>thr); end
end
function idx = ismember_runs(rid,rs), idx=[]; for r=rs(:)', idx=[idx; find(rid==r)]; end, end
function q = quantile_(x,p), x=sort(x(:)); n=numel(x); q=x(min(n,max(1,round(p*n)))); end
function y = winmean(x,rid,W)
  y = nan(size(x));
  for r = unique(rid)'
    ii = find(rid==r); v = x(ii);
    for n = W:numel(v), y(ii(n)) = mean(v(n-W+1:n)); end
  end
end
cells = {'g10','alt500','alt1000','alt2000','alt3500','alt5000','alt7000','alt10000','loit_circle','loit_los2000','loit_los5000','loit_los10000'};
small=[true true true false false false]; Nboot=500; W=12; rng(0);
fid=fopen('sweep5.csv','w'); fprintf(fid,'cell,k,pd,lo,hi,auc,E_h0,E_atk,ratio_h0,ratio_atk\n');
for c = cells
  L=load(sprintf('cell_%s.mat',c{1})); Pfa=L.Pfa; Nruns=L.Nruns;
  assert(L.done_runs==Nruns && L.done_scen==2, ['incomplete ' c{1}]);
  T0=L.T{1,1}; T1=L.T{1,2};
  if size(T0,2)<5, T0(:,5)=abs(T0(:,3)-1); T1(:,5)=abs(T1(:,3)-1); end
  T0(:,6)=winmean(T0(:,4),L.RID{1,1},W); T1(:,6)=winmean(T1(:,4),L.RID{1,2},W);
  line = sprintf('%-14s |', c{1});
  for k=1:6
    pd=pd_at_pfa(T0(:,k),T1(:,k),Pfa,small(k)); pb=zeros(Nboot,1);
    for b=1:Nboot
      rs=randi(Nruns,Nruns,1);
      pb(b)=pd_at_pfa(T0(ismember_runs(L.RID{1,1},rs),k),T1(ismember_runs(L.RID{1,2},rs),k),Pfa,small(k));
    end
    ci=quantile_(pb,[0.025 0.975]);
    t0=T0(:,k); t1=T1(:,k); t0=t0(~isnan(t0)); t1=t1(~isnan(t1)); if small(k), t0=-t0; t1=-t1; end
    [~,ord]=sort([t0;t1]); rk=zeros(size(ord)); rk(ord)=1:numel(ord);
    auc=(sum(rk(numel(t0)+1:end))-numel(t1)*(numel(t1)+1)/2)/(numel(t0)*numel(t1));
    fprintf(fid,'%s,%d,%g,%g,%g,%g,%g,%g,%g,%g\n',c{1},k,pd,ci(1),ci(2),auc,median(L.E{1,1}),median(L.E{1,2}),median(T0(:,3)),median(T1(:,3)));
    line=[line sprintf(' %.2f[%.2f,%.2f] a%.2f |',pd,ci(1),ci(2),auc)];
  end
  printf('%s E %.3f->%.3f ratio %.3f->%.3f\n', line, median(L.E{1,1}), median(L.E{1,2}), median(T0(:,3)), median(T1(:,3)));
end
fclose(fid); disp('AGGREGATED5');
