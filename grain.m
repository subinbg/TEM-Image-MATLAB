function varargout = grain(varargin)

% Last Modified by GUIDE v2.5 07-Nov-2016 02:15:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @grain_OpeningFcn, ...
                   'gui_OutputFcn',  @grain_OutputFcn, ...
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


% --- Executes just before grain is made visible.
function grain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to grain (see VARARGIN)
set(handles.text10,'String',[num2str(get(handles.slider1,'Value')) '%']);
% Choose default command line output for grain
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes grain wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = grain_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.jpg';'*.png';'*.bmp';'*.*'}, 'File Selector');
if filename == 0
    display('Canceled');
    return;
end

img = imread(fullfile(pathname,filename));
axes(handles.axes1);
imshow(img);

handles.img = img;
[handles.pathname, handles.filename, handles.ext] = fileparts(fullfile(pathname,filename));

handles.output = hObject;
guidata(hObject, handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text2.
function text2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% INPUTS
%A=imread('image2.png');
A = handles.img;
Resolution=get(handles.edit1,'String'); % micron/pixel
Resolution = str2num(Resolution);

% CALCULATIONS
%if size(A,3)==3
A=im2bw(A,graythresh(A));
%end

% Converting white space of A into value 0(false). otherwise, 1(true)
Conn=8; % Mathwork recommended connectivity for watershed transform...
[s1,s2]=size(A);
A=~bwmorph(A,'majority',10);
Resolution = Resolution/s2;

%
Poro=sum(sum(~A))/(s1*s2); % Poro = total grain size(white in original)/image size
D=-bwdist(A,'cityblock'); % D = (-1) * distance from that pixel to nearist nonzero pixel with cityblock scheme
B=medfilt2(D,[3 3]); % B = noise filtering of D; pixels in B = median value of 3by3 matrix centered by that pixel in D
B=watershed(B,Conn); % B = verifying different grain sizes; make grain boundaries obvious
% Verifying watersheds(minima) in terms of distances...
% To confirm what's going on, use
% imshow(A); figure; aa=A; aa(imdilate(B==0,ones(3,3))) = 255; figure; imshow(aa)
Pr=zeros(s1,s2);

B(:,:)=1; %removing watershed effects... (optional)
for I=1:s1
    for J=1:s2
        if A(I,J)==0 && B(I,J)~=0 % if pixel belongs to white region of A and does not belong to grain boundaries...
            Pr(I,J)=1;
        end
    end
end
Pr=bwareaopen(Pr,9,Conn);  % removing all connected components smaller than "9" pixels with connectivity used in watershed
[Pr_L,Pr_n]=bwlabel(Pr,Conn);
%bwlabel finds regions (connected pixels) with sizes larger than Conn and
%then gives the distinct regions different numbers.
%we can "select" the regions by "Z = Pr_n == n, where 1 <= n <= Pr_n.

V=zeros(Pr_n,1);
for I=1:s1
    for J=1:s2
        if Pr_L(I,J)~=0
            V(Pr_L(I,J))=V(Pr_L(I,J))+1; %storing number of pixels of each grain in column vector V.
        end
    end
end
%disp(sum(V)/numel(handles.img));
%disp(size(handles.img));
%disp(numel(handles.img));
%disp(V);
%R=Resolution.*(V./pi).^.5; % Pore radius
R = Resolution^2 *V;
%R = sort(R);

%Outputs
handles.avg=mean(R);
addr = find(R<handles.avg*get(handles.slider1,'Value')*0.01);
for I=1:numel(addr)
    Pr_L(Pr_L==addr(I)) = 0;
end
R = R(R>=handles.avg*get(handles.slider1,'Value')*0.01);
R = sort(R);

handles.avg=mean(R);
handles.std=std(R);
set(handles.text5,'String',['Avg Size = ' num2str(handles.avg)])
set(handles.text6,'String',['Std Size = ' num2str(handles.std)])
set(handles.text8,'String',['Fraction = ' num2str(sum(R)/(numel(handles.img)*Resolution^2))])

axes(handles.axes2)
%figure('units','normalized','outerposition',[0 0 1 1])
RGB=label2rgb(Pr_L,'jet', 'w', 'shuffle');
imshow(RGB)
outputfile = [handles.filename '_output.png'];
if exist([handles.filename '_outputdata.xls'],'file')==2
    delete([handles.filename '_outputdata.xls']);
end
xlswrite([handles.filename '_outputdata'],R)
imwrite(RGB,outputfile)
handles.output = hObject;
guidata(hObject, handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text10,'String',[num2str(get(handles.slider1,'Value')) '%']);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.output = hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
