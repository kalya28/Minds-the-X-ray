function varargout = xray_gui(varargin)
% XRAY_GUI MATLAB code for xray_gui.fig
%      XRAY_GUI, by itself, creates a new XRAY_GUI or raises the existing
%      singleton*.
%
%      H = XRAY_GUI returns the handle to a new XRAY_GUI or the handle to
%      the existing singleton*.
%
%      XRAY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in XRAY_GUI.M with the given input arguments.
%
%      XRAY_GUI('Property','Value',...) creates a new XRAY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before xray_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to xray_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help xray_gui

% Last Modified by GUIDE v2.5 17-Mar-2017 14:13:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @xray_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @xray_gui_Outp
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before xray_gui is made visible.
function xray_gui_OpeningFcn(hObject, eventdata, handleut args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to xray_gui (see VARARGIN)

% Choose default command line output for xray_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes xray_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = xray_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line= handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = imgetfile();
global img;
 ImgBlurSigma = 2; 
 MinHoughPeakDistance = 5;
 HoughConvolutionLength = 40; 
 HoughConvolutionDilate = 2; 
 BreakLineTolerance = 0.25; 
 breakPointDilate = 6; 

img=imread(filename);
%img{i} = imread(list(i).name);
 if pathname
     msgbox(sprintf('Error'),'Error','Error');
     return
 end
  axes(handles.axes1);le('original image');

% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img;
 ImgBlurSigma = 2; 
img = (rgb2gray(img)); 
img = imfilter(img, fspecial('gaussian', 10, ImgBlurSigma), 'symmetric'); % Denoise
axes(handles.axes2);
imshow(img),title('filtered output');

% --- Executes during object creation, after setting all properties.
function pushbutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img;
global edgeImg;
[h1,h]=c_means(img);
axes(handles.axes7);
imshow(h),title('Fuzzy C means');
pause(0.1);
boneEdges = edge(img, 'canny');
axes(handles.axes3);
imshow(boneEdges),title('Edge detection output');

boneEdges = bwmorph(boneEdges, 'close');
edgeRegs = regionprops(boneEdges, 'Area', 'PixelIdxList');
AreaList = sort(vertcat(edgeRegs.Area), 'descend');
edgeRegs(~ismember(vertcat(edgeRegs.Area), AreaList(1:2))) = [];
edgeImg = zeros(size(img, 1), size(img,2));
edgeImg(vertcat(edgeRegs.PixelIdxList)) = 1;

% --- Executes during object creation, after setting all properties.
function pushbutton3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img;
global edgeImg;
 MinHoughPeakDistance = 5;
 HoughConvolutionLength = 40; 
 HoughConvolutionDilate = 2; 
 BreakLineTolerance = 0.25; 
 breakPointDilate = 6; 


[H,T,R] = hough(edgeImg,'RhoResolution',1,'Theta',-90:2:89.5);
maxHough = max(H, [], 1);
HoughThresh = (max(maxHough) - min(maxHough))/2 + min(maxHough);
[~, HoughPeaks] = findpeaks(maxHough,'MINPEAKHEIGHT',HoughThresh, 'MinPeakDistance', MinHoughPeakDistance);
axes(handles.axes4);
plot(T, maxHough);
hold on
plot([min(T) max(T)], [HoughThresh, HoughThresh], 'r');
plot(T(HoughPeaks), maxHough(HoughPeaks), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
hold off
xlabel('Theta Value'); ylabel('Max Hough Transform');
legend({'Max Hough Transform', 'Hough Peak Threshold', 'Detected Peak'});

if numel(HoughPeaks) > 1;
    BreakStack = zeros(size(img, 1), size(img, 2), numel(HoughPeaks));
    
    for m = 1:numel(HoughPeaks);

        boneKernel = strel('line', HoughConvolutionLength, T(HoughPeaks(m)));
        kern = double(bwmorph(boate', HoughConvolutionDilate));
        BreakStack(:,:,m) = imfilter(edgeImg, kern).*edgeImg;
    end

    
    brImg = abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0;
    [BpY, BpX] = find(abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0);
    brImg = bwmorph(brImg, 'dilate', breakPointDilate);
    brReg = regionprops(brImg, 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
        'Orientation', 'Centroid');
    brReg(vertcat(brReg.Area) ~= max(vertcat(brReg.Area))) = [];

    
    brReg.EllipseCoords = zeros(100, 2);
    t = linspace(0, 2*pi, 100);
    brReg.EllipseCoords(:,1) = brReg.Centroid(1) + brReg.MajorAxisLength/2*cos(t - brReg.Orientation);
    brReg.EllipseCoords(:,2) = brReg.Centroid(2) + brReg.MinorAxisLength/2*sin(t - brReg.Orientation);

else
    brReg = [];

end

axes(handles.axes5);
imshow(img),title('fracture detection output')
hold on
colormap('gray')
if ~isempty(brReg)
    plot(brReg.EllipseCoords(:,1), brReg.EllipseCoords(:,2), 'r');
end
hold off

h = msgbox('Fracture detected','Success');

% --- Executes during object creation, after setting all properties.
function pushbutton4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
