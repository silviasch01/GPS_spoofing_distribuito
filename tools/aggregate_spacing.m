1;
% aggregate cells with FIVE detectors (col 5 = |median(r)-1| derived from the
% saved median column when absent) + a windowed worst-pair variant (k=6,
% moving mean over W epochs within each mission); Pd + 95% CI, AUC, CSVs.
function thr = upper_quantile(x,Pfa), x=sort(x(:)); n=numel(x); thr=x(min(n,max(1,ceil((1-Pfa)*n)))); end
function thr = lower_quantile(x,Pfa), x=sort(x(:)); n=numel(x); thr=x(min(n,max(1,floor(Pfa*n)+1))); end
function pd = pd_at_pfa(t0,t1,Pfa,small)
  t0=t0(~isnan(t0)); t1=t1(~isnan(t1));
  if small, thr=lower_quantile(t0,Pfa); pd=mean(t1<thr); else, thr=upper_quantile(t0,Pfa); pd=mean(t1>thr); end
end
function idx = ismember_runs(rid,rs), idx=[]; for r=rs(:)', idx=[idx; find(rid==r)]; end, end
function q = quantile_(x,p), x=sort(x(:)); n=numel(x); q=x(min(n,max(1,round(p*n)))); end
function y = winmean(x,rid,W)   % moving mean over the last W epochs within each mission
  y = nan(size(x));
  for r = unique(rid)'
    ii = find(rid==r); v = x(ii);
    for n = W:numel(v), y(ii(n)) = mean(v(n-W+1:n)); end
  end
end

gs = str2num(getenv('GS')); if isempty(gs), gs=[5 8 10 12 15 20]; end
Nboot=500; W=12; rng(0);
names={'Jansen min','Jansen max','Jansen majority','Worst-pair (two-sided)','Two-sided majority','Worst-pair, 1-min window','Two-sided majority, 1-min window'};
small=[true true true false false false false];
fid=fopen('pd4_ci.csv','w'); fe=fopen('err4_runs.csv','w'); fa=fopen('auc4.csv','w');
for g=gs
  L=load(sprintf('cell_g%d.mat',g)); Pfa=L.Pfa; Nruns=L.Nruns;
  assert(L.done_runs==Nruns && L.done_scen==2, 'cell incomplete');
  T0=L.T{1,1}; T1=L.T{1,2};
  if size(T0,2)<5, T0(:,5)=abs(T0(:,3)-1); T1(:,5)=abs(T1(:,3)-1); end
  T0(:,6)=winmean(T0(:,4),L.RID{1,1},W); T1(:,6)=winmean(T1(:,4),L.RID{1,2},W);
  T0(:,7)=winmean(T0(:,5),L.RID{1,1},W); T1(:,7)=winmean(T1(:,5),L.RID{1,2},W);
  for k=1:7
    pd=pd_at_pfa(T0(:,k),T1(:,k),Pfa,small(k));
    pb=zeros(Nboot,1);
    for b=1:Nboot
      rs=randi(Nruns,Nruns,1);
      pb(b)=pd_at_pfa(T0(ismember_runs(L.RID{1,1},rs),k),T1(ismember_runs(L.RID{1,2},rs),k),Pfa,small(k));
    end
    ci=quantile_(pb,[0.025 0.975]);
    t0=T0(:,k); t1=T1(:,k); t0=t0(~isnan(t0)); t1=t1(~isnan(t1));
    if small(k), t0=-t0; t1=-t1; end
    [~,ord]=sort([t0;t1]); rk=zeros(size(ord)); rk(ord)=1:numel(ord);
    auc=(sum(rk(numel(t0)+1:end))-numel(t1)*(numel(t1)+1)/2)/(numel(t0)*numel(t1));
    fprintf(fid,'%d,%d,%g,%g,%g\n',k,g,pd,ci(1),ci(2)); fprintf(fa,'%d,%d,%g\n',k,g,auc);
    printf('g=%3d | %-26s | Pd %.3f  [%.3f, %.3f] | AUC %.3f\n', g, names{k}, pd, ci(1), ci(2), auc);
    v0=T0(:,k); v1=T1(:,k);
    dlmwrite(sprintf('roc4_g%d_k%d_s1.csv',g,k),v0(~isnan(v0))); dlmwrite(sprintf('roc4_g%d_k%d_s2.csv',g,k),v1(~isnan(v1)));
  end
  for r=1:Nruns
    e0=median(L.E{1,1}(L.RID{1,1}==r)); e1=median(L.E{1,2}(L.RID{1,2}==r));
    fprintf(fe,'%d,%d,%g,%g\n',g,r,e0,e1);
  end
  printf('g=%3d | distortion H0 %.2f%%  attack %.2f%% | median(r) H0 %.4f attack %.4f | theta>10 frac %.3f\n', g, 100*median(L.E{1,1}), 100*median(L.E{1,2}), median(T0(:,3)), median(T1(:,3)), mean(L.theta_sat(:)));
end
fclose(fid); fclose(fe); fclose(fa); disp('AGGREGATED4');
