function varargout = demo(varargin)
% DEMO MATLAB code for demo.fig
%      DEMO, by itself, creates a new DEMO or raises the existing
%      singleton*.
%
%      H = DEMO returns the handle to a new DEMO or the handle to
%      the existing singleton*.
%
%      DEMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEMO.M with the given input arguments.
%
%      DEMO('Property','Value',...) creates a new DEMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before demo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to demo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help demo

% Last Modified by GUIDE v2.5 06-Oct-2014 12:02:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @demo_OpeningFcn, ...
                   'gui_OutputFcn',  @demo_OutputFcn, ...
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


% --- Executes just before demo is made visible.
function demo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to demo (see VARARGIN)

% Choose default command line output for demo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes demo wait for user response (see UIRESUME)
% uiwait(handles.figure1);

shnull(handles)
addpath('includes');


% --- Outputs from this function are returned to the command line.
function varargout = demo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in simButton.
function simButton_Callback(hObject, eventdata, handles)
% hObject    handle to simButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
line_bits = str2num(get(handles.linebitsEdit, 'String'));
lines_count = str2num(get(handles.linescountEdit, 'String'));
seed = str2num(get(handles.prngEdit, 'String'));

settings = default_settings();

allocEnum = get(handles.allocPopup, 'Value');
selectionEnum = get(handles.selectionPopup, 'Value');
fusionEnum = get(handles.fusionPopup, 'Value');
filtEnum = get(handles.filtlevelsPopup, 'Value');
borderEnum = get(handles.bordersPopup, 'Value');
selEnum = get(handles.selectionPopup, 'Value');

settings.bit_alloc_metric = allocEnum;
settings.disable_borders = borderEnum == 2;
settings.filter_levels = 6 - filtEnum;
settings.fusion_type = fusionEnum - 1;
settings.display = true;
if selEnum == 1
    settings.profile_mode = 0;
elseif selEnum == 2
    settings.profile_mode = -1;
elseif selEnum == 3
    settings.profile_mode = 4;
elseif selEnum == 5
    settings.profile_mode = 6;
elseif selEnum == 6
    settings.profile_mode = 8;
else
    settings.profile_mode = 0;
end

I = im2double(imread(handles.image_path));
if size(I,3) > 1
    I = rgb2gray(I);
end

oldcursor = get(handles.figure1, 'Pointer');
set(gcf, 'Pointer', 'watch');
drawnow;
[handles.imgref, ~, ~, ~, handles.imgscan] = sim_IF_reconstruction(I, line_bits, lines_count, seed, settings);
set(gcf, 'Pointer', oldcursor);

populate_available_data(handles);
show_image('Simulated reconstruction', handles);
guidata(hObject, handles);

function linebitsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to linebitsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of linebitsEdit as text
%        str2double(get(hObject,'String')) returns contents of linebitsEdit as a double
verifyEdit(hObject, 128, 16, 768, 0);
curVal = str2num(get(hObject, 'String'));
if mod(curVal, 16) ~= 0
    set(hObject, 'BackgroundColor', 'r');
    set(hObject, 'TooltipString', 'During embedding / recovery this value will be rounded to a multiple of 16.');
else
    set(hObject, 'BackgroundColor', get(0,'DefaultUicontrolBackgroundColor'));
    set(hObject, 'TooltipString', '');
end
update_capacity_info(handles)

% --- Executes during object creation, after setting all properties.
function linebitsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linebitsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function linescountEdit_Callback(hObject, eventdata, handles)
% hObject    handle to linescountEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of linescountEdit as text
%        str2double(get(hObject,'String')) returns contents of linescountEdit as a double
verifyEdit(hObject, 1024, 16, 65536, 0);

% --- Executes during object creation, after setting all properties.
function linescountEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linescountEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in allocPopup.
function allocPopup_Callback(hObject, eventdata, handles)
% hObject    handle to allocPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns allocPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from allocPopup


% --- Executes during object creation, after setting all properties.
function allocPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allocPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in selectionPopup.
function selectionPopup_Callback(hObject, eventdata, handles)
% hObject    handle to selectionPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selectionPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectionPopup


% --- Executes during object creation, after setting all properties.
function selectionPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectionPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fusionPopup.
function fusionPopup_Callback(hObject, eventdata, handles)
% hObject    handle to fusionPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fusionPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fusionPopup


% --- Executes during object creation, after setting all properties.
function fusionPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fusionPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prngEdit_Callback(hObject, eventdata, handles)
% hObject    handle to prngEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prngEdit as text
%        str2double(get(hObject,'String')) returns contents of prngEdit as a double
verifyEdit(hObject, 1234, 1, 65536, 0);

% --- Executes during object creation, after setting all properties.
function prngEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prngEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in filtlevelsPopup.
function filtlevelsPopup_Callback(hObject, eventdata, handles)
% hObject    handle to filtlevelsPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filtlevelsPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filtlevelsPopup


% --- Executes during object creation, after setting all properties.
function filtlevelsPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filtlevelsPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in bordersPopup.
function bordersPopup_Callback(hObject, eventdata, handles)
% hObject    handle to bordersPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns bordersPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from bordersPopup


% --- Executes during object creation, after setting all properties.
function bordersPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bordersPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browseButton.
function browseButton_Callback(hObject, eventdata, handles)
% hObject    handle to browseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.jpg;*.jpeg;*.png;*.tiff','Select image');
if FileName ~= 0
    set(handles.imageEditBox, 'String', FileName);
    handles.image_path = [PathName FileName];
    I = im2double(imread(handles.image_path));
    if size(I,3) > 1
        I = rgb2gray(I);
    end
    handles.img = I;
    
    set(handles.imageStatsText, 'String', sprintf('Size: %.1f Mpx, Blocks: %d', numel(I)/1024/1024, prod(size(I)/32)));
    
    set(handles.encodeButton, 'Enable', 'on');
    set(handles.simButton, 'Enable', 'on');
    populate_available_data(handles);
    show_image('Original image', handles);
    guidata(hObject, handles);
end

function imageEditBox_Callback(hObject, eventdata, handles)
% hObject    handle to imageEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageEditBox as text
%        str2double(get(hObject,'String')) returns contents of imageEditBox as a double


% --- Executes during object creation, after setting all properties.
function imageEditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in viewPopup.
function viewPopup_Callback(hObject, eventdata, handles)
% hObject    handle to viewPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns viewPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from viewPopup
list=get(handles.viewPopup,'String');
if ischar(list)
    str = list;
else
    str = list{get(handles.viewPopup,'Value')};
end

comment = '';
if all(strcmp(str,'Original image'))
    imgt = handles.img;
end

if all(strcmp(str, 'Scanned Profiles'))
    imgt = handles.imgscan;
end

if all(strcmp(str,'Simulated reconstruction'))
    imgt = handles.imgref;
end

if all(strcmp(str,'Watermarked image'))
    imgt = handles.imgw;
end

if all(strcmp(str,'Tampering source'))
    imgt = handles.imgts;
end

if all(strcmp(str,'Tampered image'))
    imgt = handles.imgt;
end

if all(strcmp(str,'Restored image'))
    imgt = handles.imgr;
end

if all(strcmp(str,'Tampering map'))
    imgt = handles.img_tampering;
end

if exist('imgt', 'var')

    if isfield(handles, 'img')
        if numel(handles.img) == numel(imgt)
            set(handles.scoresText, 'String',sprintf('%s : PSNR : %.2f dB, SSIM: %.2f %s', str, calc_psnr(im2uint8(handles.img), im2uint8(imgt)), ssim_index(im2uint8(handles.img), im2uint8(imgt)), comment));
        else
            set(handles.scoresText, 'String',sprintf('%s : Scores not available', str));
        end
    else
        set(handles.scoresText, 'String', sprintf('%s : Scores not available', str));
    end    
    
    axes(handles.axes1);
    if max(imgt(:)) > 1
        imgt = im2double(imgt);
    end
    imagesc(imgt, [0 1]); set(gca,'YTick',[]); set(gca,'XTick',[]);
    colormap gray;
    
    if get(handles.gridBox, 'Value') == 1        
        hold on;
        for x = 32:32:size(imgt,2)-1
            plot([x x], [0 size(imgt, 1)], 'y:')
        end
        for y = 32:32:size(imgt,1)-1
            plot([0 size(imgt, 2)], [y y], 'y:')
        end
        hold off;
    end

else
    axes(handles.axes1);
    shnull(handles);
end


% --- Executes during object creation, after setting all properties.
function viewPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viewPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function populate_available_data(handles)
str = cell(1,1);
found = 0;
found = found + 1;
str{found} = '<none>';
if isfield(handles, 'img')
    found = found + 1;
    str{found} = 'Original image';
end
if isfield(handles, 'imgscan')
    found = found + 1;
    str{found} = 'Scanned Profiles';
end
if isfield(handles, 'imgref')
    found = found + 1;
    str{found} = 'Simulated reconstruction';
end
if isfield(handles, 'imgw')
    found = found + 1;
    str{found} = 'Watermarked image';
end
if isfield(handles, 'imgts')
    found = found + 1;
    str{found} = 'Tampering source';
end
if isfield(handles, 'imgt')
    found = found + 1;
    str{found} = 'Tampered image';
end
if isfield(handles, 'img_tampering')
    found = found + 1;
    str{found} = 'Tampering map';
end
if isfield(handles, 'imgr')
    found = found + 1;
    str{found} = 'Restored image';
end
set(handles.viewPopup,'String',str);

function show_image(image, handles)

list = get(handles.viewPopup,'String'); 
for i = 1:length(list)
    if all(strcmp(list{i},image))
        set(handles.viewPopup, 'Value', i);
        viewPopup_Callback(handles.viewPopup, 0, handles);
    end
end

function shnull(handles)

img_null = imresize([0.5 0.6; 0.6 0.5], [16 16], 'nearest');
img_null = repmat(img_null, 32, 32);
colormap gray;
imagesc(img_null, [0,1]); set(gca,'YTick',[]); set(gca,'XTick',[]);
set(handles.scoresText, 'String', '');

function verifyEdit(hObject, default, minval, maxval, allow_float)
curVal = str2num(get(hObject, 'String'));
if isempty(curVal) 
    set(hObject, 'String', num2str(default));
    return
end
if curVal < minval
    set(hObject, 'String', num2str(minval));
    curVal = minval;
end
if curVal > maxval
    set(hObject, 'String', num2str(maxval));
    curVal = maxval;
end
if allow_float == false && curVal ~= floor(curVal)
    curVal = floor(curVal);
    set(hObject, 'String', num2str(curVal));
end

function update_capacity_info(handles)
curVal = str2num(get(handles.linebitsEdit, 'String'));
curVal = round(curVal / 16);
set(handles.capacityText, 'String', sprintf('Cap.: %d bpb (%.3f bpp)', curVal, curVal/64))
firstCoeff = str2num(get(handles.firstCoeffEdit, 'String'));
set(handles.lastCoeffEdit, 'String', num2str(firstCoeff + curVal - 1));

% --- Executes on button press in gridBox.
function gridBox_Callback(hObject, eventdata, handles)
% hObject    handle to gridBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gridBox
viewPopup_Callback(handles.viewPopup, 0, handles);



function quantStepEdit_Callback(hObject, eventdata, handles)
% hObject    handle to quantStepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of quantStepEdit as text
%        str2double(get(hObject,'String')) returns contents of quantStepEdit as a double
verifyEdit(hObject, 0.055, 0.02, 0.15, 1);

% --- Executes during object creation, after setting all properties.
function quantStepEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quantStepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function firstCoeffEdit_Callback(hObject, eventdata, handles)
% hObject    handle to firstCoeffEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of firstCoeffEdit as text
%        str2double(get(hObject,'String')) returns contents of firstCoeffEdit as a double
verifyEdit(hObject, 7, 2, 63, 0);
update_capacity_info(handles)

% --- Executes during object creation, after setting all properties.
function firstCoeffEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to firstCoeffEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lastCoeffEdit_Callback(hObject, eventdata, handles)
% hObject    handle to lastCoeffEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lastCoeffEdit as text
%        str2double(get(hObject,'String')) returns contents of lastCoeffEdit as a double


% --- Executes during object creation, after setting all properties.
function lastCoeffEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lastCoeffEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in encodeButton.
function encodeButton_Callback(hObject, eventdata, handles)
% hObject    handle to encodeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
line_bits = str2num(get(handles.linebitsEdit, 'String'));
quantization_step = str2num(get(handles.quantStepEdit, 'String'));
seed = str2num(get(handles.prngEdit, 'String'));

settings = default_settings();

allocEnum = get(handles.allocPopup, 'Value');
selectionEnum = get(handles.selectionPopup, 'Value');
fusionEnum = get(handles.fusionPopup, 'Value');
filtEnum = get(handles.filtlevelsPopup, 'Value');
borderEnum = get(handles.bordersPopup, 'Value');

settings.bit_alloc_metric = allocEnum;
settings.disable_borders = borderEnum == 2;
settings.filter_levels = 6 - filtEnum;
settings.fusion_type = fusionEnum - 1;
settings.display = 1;

I = im2double(handles.img);
if size(I,3) > 1
    I = rgb2gray(I);
end

oldcursor = get(handles.figure1, 'Pointer');
set(gcf, 'Pointer', 'watch');
drawnow;
handles.imgw = encode(I, line_bits, quantization_step, seed, settings);
set(gcf, 'Pointer', oldcursor);

set(handles.saveProtButon, 'Enable', 'on');
populate_available_data(handles);
show_image('Watermarked image', handles);
guidata(hObject, handles);

% --- Executes on button press in decodeButton.
function decodeButton_Callback(hObject, eventdata, handles)
% hObject    handle to decodeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
line_bits = str2num(get(handles.linebitsEdit, 'String'));
quantization_step = str2num(get(handles.quantStepEdit, 'String'));
seed = str2num(get(handles.prngEdit, 'String'));

settings = default_settings();

allocEnum = get(handles.allocPopup, 'Value');
selectionEnum = get(handles.selectionPopup, 'Value');
fusionEnum = get(handles.fusionPopup, 'Value');
filtEnum = get(handles.filtlevelsPopup, 'Value');
borderEnum = get(handles.bordersPopup, 'Value');

settings.bit_alloc_metric = allocEnum;
settings.disable_borders = borderEnum == 2;
settings.filter_levels = 6 - filtEnum;
settings.fusion_type = fusionEnum - 1;
settings.display = 1;

I = im2double(handles.imgt);
if size(I,3) > 1
    I = rgb2gray(I);
end

oldcursor = get(handles.figure1, 'Pointer');
set(gcf, 'Pointer', 'watch');
drawnow;
[handles.imgr, t_map] = decode(I, line_bits, seed, settings);
set(gcf, 'Pointer', oldcursor);

handles.img_tampering = overlay_tampering(I, t_map);

set(handles.statsValidText, 'String', sprintf('Valid lines: %d (%.1f%%)', sum(t_map(:) == 0), 100*sum(t_map(:) == 0)/numel(t_map)));
set(handles.statsCrcText, 'String', sprintf('Damaged CRC: %d (%.1f%%)', sum(t_map(:) == 1), 100*sum(t_map(:) == 1)/numel(t_map)));
set(handles.statsSvmText, 'String', sprintf('SVM rejected: %d (%.1f%%)', sum(t_map(:) == 2), 100*sum(t_map(:) == 2)/numel(t_map)));

set(handles.saveAuthButton, 'Enable', 'on');
populate_available_data(handles);
show_image('Restored image', handles);
guidata(hObject, handles);

% --- Executes on button press in saveProtButon.
function saveProtButon_Callback(hObject, eventdata, handles)
% hObject    handle to saveProtButon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile('*.png;','Save protected image as...');
if FileName ~= 0
    imwrite(handles.imgw, [PathName FileName]);
end

% --- Executes on button press in saveAuthButton.
function saveAuthButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveAuthButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile('*.png;','Save restored image as...');
if FileName ~= 0
    imwrite(handles.imgr, [PathName FileName]);
end



function imageInspectedEdit_Callback(hObject, eventdata, handles)
% hObject    handle to imageInspectedEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageInspectedEdit as text
%        str2double(get(hObject,'String')) returns contents of imageInspectedEdit as a double


% --- Executes during object creation, after setting all properties.
function imageInspectedEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageInspectedEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in inspectedBrowseButton.
function inspectedBrowseButton_Callback(hObject, eventdata, handles)
% hObject    handle to inspectedBrowseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.jpg;*.jpeg;*.png;*.tiff','Select image for verification');
if FileName ~= 0
    set(handles.imageInspectedEdit, 'String', FileName);
    handles.image_inspected_path = [PathName FileName];
    I = imread(handles.image_inspected_path);
    if size(I,3) > 1
        I = rgb2gray(I);
    end
    handles.imgt = I;
    
    set(handles.imageInspStatsText, 'String', sprintf('Size: %.1f Mpx, Blocks: %d', numel(I)/1024/1024, prod(size(I)/32)));
    
    set(handles.decodeButton, 'Enable', 'on');
    populate_available_data(handles);
    show_image('Tampered image', handles);
    guidata(hObject, handles);
end


% --- Executes on button press in aboutButton.
function aboutButton_Callback(hObject, eventdata, handles)
% hObject    handle to aboutButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
about_text = {'Iterative Filtering for Semi-Fragile Self-Recovery - Interactive Demo','',...
    'Authors: Pawel Korus, Jaroslaw Bialas, Andrzej Dziech 2014', '(C) AGH University of Science and Technology 2014','',...
    'E-mail: pkorus@agh.edu.pl',...
    'Web: http://kt.agh.edu.pl/~korus','','This demonstration has been provided as supplementary materials to paper',...
    '"Iterative Filtering for Semi-Fragile Self-Recovery" presented during','IEEE International Workshop on Information Forensics and Security 2014.', '', ...
    'This demo aims to facilitate reproduction of the results from the paper.' ...
    'Personal use for educational and research purposes is permitted. Any other use, including commercial purposes, is strictly prohibited. The code is provided "as-is", without warranties of any kind.','',...
    'The research leading to these results received funding from the European ',...
    'Regional Development Fund under INSIGMA project no. POIG.01.01.02-00-062/09.'};
msgbox(about_text, 'About', 'help');
