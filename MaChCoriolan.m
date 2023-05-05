% Open the play text, read it, get the length of it.
f = fopen('coriolan.txt','r');
X = fread(f,inf,'uint8');
txt = char(X)';
fclose(f);
L = numel(txt);
% Make k-blocks of the play, which will be our states, find all unique
% states.
k = 2;
A = zeros(L-k,k,'uint8');
for j=1:L-k+1
    A(j,:) = X(j:j+k-1);
end
[C,IA,IC] = unique(A,'rows');
eos = IC(end);
% Find all states with multiple unique transitions. These are the only
% states that need to be arithmetically encoded.
mon = [];
for i=1:numel(IA)
    icf = find(IC==i);
    if numel(icf) > 1
        co = unique(IC(min(icf+1,numel(IC))));
        if numel(co) > 1
            mon(numel(mon)+1) = i;
        end
    end
end

% Our emissions are the end letter of the block.
[HC, HIA, HIC] = unique(X); % All of these are possible outputs of our Markov Chain.
N1 = size(C,1); S1 = sparse(N1,N1)+1;E1 = sparse(numel(IA),numel(HC))+1; % Markov States, Making the Transition Matrix.
emlist = HIC(2:numel(HIC)-1);
for i=1:L-k
    S1(IC(i),IC(i+1)) = S1(IC(i),IC(i+1)) + 100; % Large weight to transitions that are present.
    linu = emlist(i);
    E1(IC(i),linu) = E1(IC(i),linu) + 100;
end

esttr = S1./(eps+sum(S1,2)); % Estimate of the transition probabilities.
estemit = E1./sum(E1,2);

y = {}; % States we hit.
sec = zeros(1,numel(mon)); % State Entropy Matrix.
for i=1:numel(mon)
    siq = mon(i); % State to encode.
    icf = find(IC==siq);
    if siq==eos
        [Cs, IAs, ICs] =  unique(IC(icf(1:numel(icf)-1)+1));
    else
        [Cs, IAs, ICs] =  unique(IC(icf+1));
    end
    % Find possible states.    
    sq = S1(siq,Cs);sq2 = esttr(siq,Cs); % States we hit, probabilities of states.
    sec(i) = -sum(sq2.*log2(sq2)); % Entropy of the state we are in.
    y{i} = Cs; % State list.
    code{i} = arithenco(ICs,sq); % Code for state list.
end

% Finding the stationart probability of states.
[V,D] = eigs(esttr',6);
D=diag(D);
Da = abs(D);

[~,idx] = max(Da);
mu = V(:,idx)';                           
mu = mu./sum(mu);

% Theoretical Entropy
muentst = mu(mon);
H = muentst*sec';

cls = 0;
for i=1:numel(code)
    cls = cls + numel(code{i});
end
% Arithmetical Code Bits Per Symbol
acb = cls/numel(IC);
%% HMM Coding

q = L/10;r = randi([1, q], 1, 400);
seqfeed = {};
for i=1:400
    seqfeed{i,1} = emlist(10*(r(i)-1):10*r(i)-1)'; 
end
[trref, emitref] = hmmtrain(seqfeed, esttr, estemit, "verbose", true,"MaxIterations",5);
% The emission probabilities imply that each state emits is strongly
% predicted to output one output only. So only encode the state transitions.
for i=1:numel(mon)
    siq = mon(i);sq = ceil(10000*trref(siq,:))+1;sq2 = trref(siq,:)+0.0001;sq2 = sq2/sum(sq2);
    icf = find(IC==siq);
    if siq==eos
        Cs = IC(icf(1:numel(icf)-1)+1);
    else
        Cs = IC(icf+1);
    end
    sec(i) = -sum(sq2.*log2(sq2)); % Entropy of the state we are in.
    y2{i} = Cs;
    code2{i} = arithenco(Cs,sq);
end

% Theoretical Entropy
[V,D] = eigs(trref',6);
D=diag(D);
Da = abs(D);

[~,idx] = max(Da);
mu = V(:,idx)';                           
mu = mu./sum(mu);

muentst2 = mu(mon);
H2 = muentst2*sec';

cls2 = 0;
for i=1:numel(code)
    cls2 = cls2 + numel(code2{i});
end
% Arithmetical Code Bits Per Symbol
acb2 = cls2/numel(IC);
% Weirdly this is smaller than the theoretical, but the theoretical is not perfectly calculated.