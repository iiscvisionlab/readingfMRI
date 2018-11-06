%manyimagesc       -> view many images, contained in a cell array 
% h = manyimagesc(image_array,titles,nx,ny)
% Required inputs
%    image_array   = cell array or a 3d matrix containing images
% Optional inputs
%    titles        = cell/numeric array containing titles for each image
%    nx            = number of subplot rows
%    ny            = number of subplot columns 

%  SP Arun
%  Changelog
%      5/1/2007 (SPA) First version

function h = manyimagesc(image_array,titles,nx,ny,whiteflag,figh)
if(~iscell(image_array)) % then treat as 3-d matrix
    for i = 1:size(image_array,3)
        A{i} = image_array(:,:,i); 
    end
    image_array = A; 
end
if(nargin<3)
    nx = ceil(sqrt(length(image_array))); ny = nx; 
end
if(~exist('whiteflag')),whiteflag = 0; end; 

titleflag=0; % by default no titles
if(nargin>1)
    if(~exist('titles')|isempty(titles))
        titles = '';
    else
        titleflag = 1; 
    end
end

if(exist('figh')), h(1) = figh; else,h(1)=figure; end
for i = 1:length(image_array)
    h(i+1)=subplot(nx,ny,i);
    if(whiteflag==0)
        imagesc(image_array{i}); 
    else
        imagesc(255-image_array{i}); 
    end
    axis image; colormap gray; axis off;
    set(gca,'XTick',[],'YTick',[],'FontSize',8);
    if(titleflag==1 & isempty(titles))
        title(h(i+1),sprintf('image %d',i));
    elseif(titleflag==1 & iscell(titles))
        title(h(i+1),titles{i});
    elseif(titleflag==1 & isnumeric(titles))
        title(h(i+1),sprintf('%3.3g',titles(i)));
    end
end

return