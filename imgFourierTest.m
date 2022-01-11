clc; clear;

img = imread('img1.png','BackgroundColor',[0 0 0]);
pxl = img(:,:,1) ~= 255;

[r,c] = find(pxl==1,1);
f(1,:) = [r,c];
pxl(r,c) = 0;
i = -1:1;

j = 2;
while 1
    near = pxl(r+i,c+i);
    near(2,2) = 0;
    [rt,ct] = find(near==1,1);
    if isempty(rt)
        break;
    end
    r = r+i(rt);
    c = c+i(ct);
    f(j,:) = [r,c];
    pxl(r,c) = 0;
    j = j+1;
end

Mid = size(pxl)/2;

f = f-Mid;

t = linspace(0,1,size(f,1))';

N = 100;
cp = zeros(1,N/2);
cn = zeros(1,N/2);
c0 = sum(f(:,1)+1i*f(:,2))/size(f,1);
for k = 1:N/2
    cp(k) = sum((f(:,1)+1i*f(:,2)).*exp(-k*2i*pi*t))/size(f,1); %positive coefficients
    cn(k) = sum((f(:,1)+1i*f(:,2)).*exp(k*2i*pi*t))/size(f,1); %negative coefficients
end

v = VideoWriter('test.avi');
open(v);

figure(1); clf;
plot(f(:,1),f(:,2))
hold on;
quiver(0,0,real(c0),imag(c0),'AutoScale','off')

for tt = t'
    vecSum = c0;

    if tt == 0
        for k = 1:N/2
            quiver(real(vecSum),imag(vecSum),real(cp(k)*exp(k*2i*pi*tt)),imag(cp(k)*exp(k*2i*pi*tt)),'AutoScale','off')
            vecSum = vecSum+cp(k)*exp(k*2i*pi*tt);
            quiver(real(vecSum),imag(vecSum),real(cn(k)*exp(-k*2i*pi*tt)),imag(cn(k)*exp(-k*2i*pi*tt)),'AutoScale','off')
            vecSum = vecSum+cn(k)*exp(-k*2i*pi*tt);
        end
    else
        a = gca;
        for k = 1:N/2
            a.Children(end-2*k).XData = real(vecSum);
            a.Children(end-2*k).YData = imag(vecSum);
            a.Children(end-2*k).UData = real(cp(k)*exp(k*2i*pi*tt));
            a.Children(end-2*k).VData = imag(cp(k)*exp(k*2i*pi*tt));
            vecSum = vecSum+cp(k)*exp(k*2i*pi*tt);
            a.Children(end-2*k-1).XData = real(vecSum);
            a.Children(end-2*k-1).YData = imag(vecSum);
            a.Children(end-2*k-1).UData = real(cn(k)*exp(-k*2i*pi*tt));
            a.Children(end-2*k-1).VData = imag(cn(k)*exp(-k*2i*pi*tt));
            vecSum = vecSum+cn(k)*exp(-k*2i*pi*tt);
        end
    end
    frame = getframe(gcf);
    writeVideo(v,frame);
    %drawnow;
end
close(v)
