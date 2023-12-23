function w = calclatwts(nlat, sig, rad)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

%nlat is the size of the window. It must be an odd number.
%sig: max strength of the connections
%rad: radius of neighborhood

%% CODE
w = zeros(nlat, nlat);
ic = (nlat+1)./2;
jc = (nlat+1)./2;


for i = 1:nlat,
   for j = 1:nlat,
        dis = (i-ic).*(i-ic) + (j-jc).*(j-jc);
        w(i, j) = sig.*exp(-dis./(rad.*rad)) ; 
        if(i==j) %&& (i==ic) && (i==jc)
            w(i,j)=0;
        end
   end
end

% w(round(nlat/2),round(nlat/2))=0;

