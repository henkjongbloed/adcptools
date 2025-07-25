function [ADCPout]=readADCP(files,varargin)
%READADCP Reads a RDI adcp binary data file
%         [ADCP]=readADCP(FILENAME) reads adcp data and outputs all data
%         
%         FILENAME can be a character array with filenames on each row or a
%         cell of strings containing filenames in each element. The
%         function will read multiple files even when configuration differs
%         amongst them. In this case the function will issue a warning.
%         Settings will be stored once for each file
%
%         [ADCP]=readADCP(FILENAME,VARS) reads adcp data given in VARS
%         VARS should be a character array containing one of the following
%         characters:
%         e   ensemble info
%         v   velocity data
%         c   correlation data
%         h   echo intensity
%         p   percentage good
%         b   bottom track data
%         x   external data stored in the adcp file
%         f   external data stored in external winriver ascii-files (The
%             function will look for the files based on the raw-data
%             filenames. Being implemented!
%
%         Example: ADCP=readADCP('data_000.000; data_001.000','vEyx') will
%         read velocity, ensemble information and external data from the
%         two files

%    Copyright 2007-2021 Bart Vermeulen, David Vermaas, Maximiliano Sassi,
%       Judith Zomer
%
%    This file is part of ADCPTools.
%
%    ADCPTools is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    ADCPTools is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with ADCPTools.  If not, see <http://www.gnu.org/licenses/>.


global dataout WRII_GGA_ID WRII_HDT_ID WRII_VTG_ID WRII_DBT_ID             %Define global container for all the data

WRII_GGA_ID=104;
WRII_HDT_ID=107;
WRII_VTG_ID=101;
WRII_DBT_ID=102;



%% Check input
if verLessThan('matlab','7.13')
    error(nargchk(1,2,nargin)) %#ok<NCHKN>
    error(nargoutchk(0,1,nargout)) %#ok<NCHKE>
else
    narginchk(1,2)
    nargoutchk(0,1)
end
inp=inputParser;                                                           % Create an object of the InputParser class
inp.addRequired('files',@(x) (isstring(x) | iscellstr(x) | ischar(x)));                  % Add the required variable 'files' and check for right format
inp.addOptional('flags','vhcpbexa',@ischar);                                % Add the optional argument 'flags' and check for right format
inp.parse(files,varargin{:});                                              % Parse input

if ischar(files)
    files=cellstr(files);                                                % Change character array to cell
end
flags=inp.Results.flags;
clear inp

%Check filenames
nfiles=numel(files);
enscnt=0;
nValidFiles=0;                                                             % Set valid files counter to 0
fileid=[];
for cntfile=1:nfiles                                                       % Loop for files
    [fid,openmessage]=fopen(files{cntfile},'r');                           % Open file in read-only binary mode
    if fid==-1                                                             % If opening file fails
        warning('readADCP:WrongFile',[openmessage,': ',files{cntfile}])   % Show warning
        continue                                                           % Continue to next file
    end

%% Find data position-information in file
    fpos=0;                                                                %Set file position pointer to beginning of file
    while 1                                                                %Loop
        headpos=searchHead(fid,fpos);                                      %Search for header
        if headpos==-1, break, end                                         %If no header was found, exit the loop
        if ~isvalens(fid,headpos)                                          %If ensemble is not valid
           warning('readADCP:CorruptEnsemble',['Invalid ensemble: ',num2str(enscnt+1)]) %Generate warning
           fpos=headpos+1;                                                 %Point fpos to position after last header
           continue                                                        %Restart loop
        end
        fseek(fid,headpos+2,-1);                                           %Move to position to read ensemble size
        EnsBytes=fread(fid,1,'uint16=>double');                            %read number of bytes in ensemble
        fpos=headpos+EnsBytes+2;                                           %Point to byte behind last read ensemble
        enscnt=enscnt+1;                                                   %Increase ensemble counter with one
        fileid(enscnt)=fid;%#ok<AGROW>                                     %Store fileid for each ensemble  
        EnsStart{enscnt}=headpos;%#ok<AGROW>                               %Store starting position of ensemble %#ok<AGROW> 
        [~,DataOffset{enscnt},DataHeader{enscnt}]=...
            readhead(fid,headpos);%#ok<AGROW>                              %Read info about where to find which data %#ok<AGROW> 
    end
    if ~any(fileid==fid)                                                   %If no valid ensembles are found in current file
        warning('readADCP:NoEnsembleInFile',...
            ['No valid ensembles in file: ',files{cntfile}])               %Generate a warning
        fclose(fid);                                                       %Close file
    else
        disp(['Found ',num2str(length(find(fileid==fid))),...
            ' valid ensembles in file: ',files{cntfile}]);                 %Display amount of valid ensembles in current file
        nValidFiles=nValidFiles+1;                                         %Record number of valid files
        ValidFilesId(nValidFiles)=fid;%#ok<AGROW>                          %Record Id of valid files
    end
end
nens=enscnt;                                                             % Calculate total number of ensembles
if nens<1                                                                  % If no ensmble is found
    error('readADCP:NoEnsemble','Could not find any valid ensemble')      % Generate error
end
disp(['Found ',num2str(nens),' valid ensembles'])                           % Display total number of valid ensembles
clear enscnt
dataout.FileNumber=fileid;
for cntfiles=1:nValidFiles
    dataout.FileNumber(fileid==ValidFilesId(cntfiles))=cntfiles;           %Record in output which ensembles belong to which FL setting field
end

disp('Reading data...')

%% Read the fixed leader
%FixedLeader is assumed to be the same for the whole data file and is read
%only for the first ensemble. This is usually true as it will only change 
%when new commands are sent to the ADCP. This can only be done by 
%interrupting the pinging. If the Fixed leaders of the files differ a
%warning is generated. Take care concatenating files with different
%settings!!! All setting will be recorded once for each file

f = waitbar(0,'Reading data...');
dataout = initFL(nens, dataout);
for cntens = 1:nens
    Ndatablock = find(DataHeader{cntens}(:,1) == 0, 1);
    if isempty(Ndatablock)
        continue
    else
        fpos = EnsStart{cntens} + DataOffset{cntens}(Ndatablock);
        dataout = readFL(fileid(cntens), fpos, cntens, dataout);
    end
    %disp(cntens)
    waitbar(cntens/nens,f,'Reading data...');
end
close(f)

% for cntfiles=1:nValidFiles
%     dataout.FileNumber(fileid==ValidFilesId(cntfiles))=cntfiles;           %Record in output which ensembles belong to which FL setting field
%     firstens=find(fileid==ValidFilesId(cntfiles),1);                       %Determine which ensemble is the first for the current file
%     posFL=EnsStart{firstens}+...
%         DataOffset{firstens}(DataHeader{firstens}(:,1)==0,1);              %Locate the fixed leader
%     readFL(ValidFilesId(cntfiles),posFL,cntfiles);                         %Call function to read the fixed leader
% end
% if ~CheckFL(nValidFiles)                                                   % If configuration in fixed leader changes
%     warning('readADCP:ConfigChange',...                                   % Generate warning
%         'The configuration seems to change between files. \n Using readADCP on files with different configuration is not recommended, but will still work') 
% end
%% Preallocate and read data


nbins=max(dataout.nbins);                                                  %Find largest amount of bins in order to make matrices in which all data fits
AllHeaders=vertcat(DataHeader{:});                                         %Put all data headers in one vector
AllHeaders(AllHeaders(:,1)~=8226,2)=65535;
%Read velocity data
if any(regexpi(flags,'v'))                                                 % Check if velocity data must be read (that is when only 1 argument or v in second argument)
    if ~any(AllHeaders(:,1)==256)                                          % Check if any velocity data is present
        warning('readADCP:NoVelData','No velocity data found.')           % Warn if no velocity data is present
    else                                                                   % Otherwise
        dataout.VEL=ones(nbins,nens,4,'int16')*-32768;                     % Preallocate velocity data (null value by default)
        for cntens=1:nens                                                  % Loop for all ensembles
            CurrentNbins=dataout.nbins(dataout.FileNumber(cntens));        % Read the number of bins in current file
            Ndatablock=find(DataHeader{cntens}(:,1)==256, 1);              % Search for velocity data in current ensemble
            if isempty(Ndatablock)                                         % If data is not found
                continue                                                   % Leave matrix empty (actually full of zeros, change to -32768???)
            else                                                           % Otherwise
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);      % Calculate position of velocity data in file
                dataout.VEL(1:CurrentNbins,cntens,:)=...
                    reshape(readVEL(fileid(cntens),...
                    fpos,CurrentNbins),[],1,4);                            % Read velocity data of current ensemble
            end                                                            % Note: for explanation on velocity data see under the reading function
        end
    end
end

%Read Echo intensity
if any(regexpi(flags,'h'))                                                 %Same as above but now for echo intensity data
    if ~any(AllHeaders(:,1)==768)
        warning('readADCP:NoEchoData','No echo intensity data found.')
    else
        dataout.ECHO=zeros(nbins,nens,4,'uint8');
        for cntens=1:nens
            CurrentNbins=dataout.nbins(dataout.FileNumber(cntens));        % Read the number of bins in current file
            Ndatablock=find(DataHeader{cntens}(:,1)==768, 1);
            if isempty(Ndatablock)
                continue
            else
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
                dataout.ECHO(1:CurrentNbins,cntens,:)=...
                    reshape(readECHO(fileid(cntens),...
                    fpos,CurrentNbins),[],1,4);
            end
        end
    end
end

% Read correlation magnitude data
if any(regexpi(flags,'c'))                                                 %Same as above but now for correlation magnitude data
    if ~any(AllHeaders(:,1)==512)
        warning('readADCP:NoCorrData','No correlation magnitude data found.')
    else
        dataout.CORR=zeros(nbins,nens,4,'uint8');
        for cntens=1:nens
            CurrentNbins=dataout.nbins(dataout.FileNumber(cntens));            % Read the number of bins in current file
            Ndatablock=find(DataHeader{cntens}(:,1)==512, 1);
            if isempty(Ndatablock)
                continue
            else
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
                dataout.CORR(1:CurrentNbins,cntens,:)=...
                    reshape(readCORR(fileid(cntens),...
                    fpos,CurrentNbins),[],1,4);
            end
        end
    end
end

%Read percentage good data
if any(regexpi(flags,'p'))                              %Same as above but now for percentage good data
    if ~any(AllHeaders(:,1)==1024)
        warning('readADCP:NoPercData','No percentage good data found.')
    else
        dataout.PERC=zeros(nbins,nens,4,'uint8');
        for cntens=1:nens
            CurrentNbins=dataout.nbins(dataout.FileNumber(cntens));            % Read the number of bins in current file
            Ndatablock=find(DataHeader{cntens}(:,1)==1024, 1);
            if isempty(Ndatablock)
                continue
            else
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
                dataout.PERC(1:CurrentNbins,cntens,:)=reshape(readPERC(fileid(cntens),fpos,CurrentNbins),[],1,4);
            end
        end
    end
end

%Read ambient sound data
if any(regexpi(flags,'a'))                              %
    if ~any(AllHeaders(:,1)==3074)                      % check if same for other data files
        warning('readADCP:NoASData','No ambient sound data found.')
    else
        dataout.ambsound=zeros(nens,4,'uint8'); 
        for cntens=1:nens 
            Ndatablock=find(DataHeader{cntens}(:,1)==3074, 1);
            if isempty(Ndatablock)
                continue
            else
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
                dataout.ambsound(cntens,:)=readAS(fileid(cntens),fpos);
            end
        end
    end
end


%Read bottom track data
if any(regexpi(flags,'b'))                              %Same as above but now for bottom track data
    if ~any(AllHeaders(:,1)==1536)
        warning('readADCP:NoBtData','No bottom track data found.')
    else
        initBT(nens);
        for cntens=1:nens
            Ndatablock=find(DataHeader{cntens}(:,1)==1536, 1);
            if isempty(Ndatablock)
                continue
            else
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
                readBT(fileid(cntens),fpos,cntens);
            end
        end
    end
end

%Read stream pro leader format
if any(AllHeaders(:,1)==20480)
    initSP(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==20480, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readSP(fileid(cntens),fpos,cntens);
        end
    end
end
%Read variable leader data
if any(regexpi(flags,'e'))                              %Same as above but now for variable leader data
    if ~any(AllHeaders(:,1)==128)
        warning('readADCP:NoVLData','No ensemble (Variable Leader) data found.')
    else
        initVL(nens);
        for cntens=1:nens
            Ndatablock=find(DataHeader{cntens}(:,1)==128, 1);
            if isempty(Ndatablock)
                continue
            else
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
                readVL(fileid(cntens),fpos,cntens);
            end
        end
    end
end

%Read instrument transformation matrix
if any(AllHeaders(:,1)==12800)
    initTM(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==12800, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readTM(fileid(cntens),fpos,cntens);
        end
    end
end

% Read vertical beam range data
if any(AllHeaders(:,1) == 16640)
    initVBeamRange(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==16640, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readVBeamRange(fileid(cntens),fpos,cntens);
        end
    end
end

% Read vertical beam profile leader data
if any(AllHeaders(:,1) == 3841)
    initVBeamLeader(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==3841, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readVBeamLeader(fileid(cntens),fpos,cntens);
        end
    end
end

% Read vertical beam velocity data
if any(AllHeaders(:,1) == 2560)
    initVBeamVelocity(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==2560, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readVBeamVelocity(fileid(cntens),fpos,cntens);
        end
    end
end

% Read vertical beam corr data
if any(AllHeaders(:,1) == 2816)
    initVBeamCorr(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==2816, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readVBeamCorr(fileid(cntens),fpos,cntens);
        end
    end
end

% Read vertical beam echo data
if any(AllHeaders(:,1) == 3072)
    initVBeamEcho(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==3072, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readVBeamEcho(fileid(cntens),fpos,cntens);
        end
    end
end

% Read vertical beam percentage good data
if any(AllHeaders(:,1) == 3328)
    initVBeamPerc(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==3328, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readVBeamPerc(fileid(cntens),fpos,cntens);
        end
    end
end

% Read vertical beam status data
if any(AllHeaders(:,1) == 3584)
    initVBeamStat(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==3584, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readVBeamStat(fileid(cntens),fpos,cntens);
        end
    end
end


% Read surface layer velocity leader
if any(AllHeaders(:,1) == 16)
    initSLayerVelLeader(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==16, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readSLayerVelLeader(fileid(cntens),fpos,cntens);
        end
    end
end

% Read surfacelayer velocity data
if any(AllHeaders(:,1) == 272)
    initSLayerVelocity(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==272, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readSLayerVelocity(fileid(cntens),fpos,cntens);
        end
    end
end

% Read surfacelayer correlation data
if any(AllHeaders(:,1) == 528)
    initSLayerCorr(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==528, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readSLayerCorr(fileid(cntens),fpos,cntens);
        end
    end
end

% Read surfacelayer echo data
if any(AllHeaders(:,1) == 784)
    initSLayerEcho(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==784, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readSLayerEcho(fileid(cntens),fpos,cntens);
        end
    end
end

% Read surfacelayer percentage  good data
if any(AllHeaders(:,1) == 1040)
    initSLayerPerc(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==1040, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readSLayerPerc(fileid(cntens),fpos,cntens);
        end
    end
end

% Read surfacelayer status data
if any(AllHeaders(:,1) == 1296)
    initSLayerStat(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==1296, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readSLayerStat(fileid(cntens),fpos,cntens);
        end
    end
end

% Read Auto Mode 3 data
if any(AllHeaders(:,1) == 17409)
    initAutoMode3(nens);
    for cntens=1:nens
        Ndatablock=find(DataHeader{cntens}(:,1)==17409, 1);
        if isempty(Ndatablock)
            continue
        else
            fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
            readAutoMode3(fileid(cntens),fpos,cntens);
        end
    end
end

% Read external data
if any(regexpi(flags,'x'))                                                 % If external data is to be read
    if any(AllHeaders(:,1)==8226)                                          % Check if any WinRiver II data is available
        if any(AllHeaders(:,2)==WRII_GGA_ID)                               % If WRII GGA data is available
            disp('WinRiver II GGA data found')                             % Display that this data is found
            dd=diff(find(AllHeaders(:,2)==WRII_GGA_ID &...
                AllHeaders(:,1)==8226));                                   % Determine distance between all WRII GGA headers
            nblocks=max(diff(find(dd~=1)));                                % Determine maximum distance (in datablocks) between non WRII GGA headers (i.e. maximum WRII GGA headers in one ensembles)
            initNMEAGGA(nens,nblocks);                                     % Initialize variables for WRII GGA data
            for cntens=1:nens                                              % Loop for every ensemble
                Ndatablock=find(DataHeader{cntens}(:,2)==WRII_GGA_ID &...
                    DataHeader{cntens}(:,1)==8226);                        % Determine Number(s) of datablock(s) with WRII GGA data
                if isempty(Ndatablock)                                     % If no WRII GGA data is available
                    continue                                               % Continue with next ensemble
                else                                                       % Otherwise
                    for cntDatBlock=1:length(Ndatablock)                   % Loop for every available datablock
                        fpos=EnsStart{cntens}+...
                            DataOffset{cntens}(Ndatablock(cntDatBlock));   % Determine position of data
                        readNMEAGGA(fileid(cntens),fpos,...
                            cntens,cntDatBlock);                           % Call function to read out the data
                    end
                end
            end
        end
        if any(AllHeaders(:,2)==WRII_HDT_ID)                               % Same as above but for WRII HDT data
            disp('WinRiver II HDT data found')
            dd=diff(find(AllHeaders(:,2)==WRII_HDT_ID...
                & AllHeaders(:,1)==8226));
            nblocks=max(diff(find(dd~=1)));
            initNMEAHDT(nens,nblocks);
            for cntens=1:nens
                Ndatablock=find(DataHeader{cntens}(:,2)==WRII_HDT_ID &...
                    DataHeader{cntens}(:,1)==8226);
                if isempty(Ndatablock)
                    continue
                else
                    for cntDatBlock=1:length(Ndatablock)
                        fpos=EnsStart{cntens}+...
                            DataOffset{cntens}(Ndatablock(cntDatBlock));
                        readNMEAHDT(fileid(cntens),fpos,cntens,cntDatBlock);
                    end
                end
            end
        end
        if any(AllHeaders(:,2)==WRII_VTG_ID)                               % Same as above but for WRII VTG data
            disp('WinRiver II VTG data found')
            dd=diff(find(AllHeaders(:,2)==WRII_VTG_ID &...
                AllHeaders(:,1)==8226));
            nblocks=max(diff(find(dd~=1)));
            initNMEAVTG(nens,nblocks);
            for cntens=1:nens
                Ndatablock=find(DataHeader{cntens}(:,2)==WRII_VTG_ID &...
                    DataHeader{cntens}(:,1)==8226);
                if isempty(Ndatablock)
                    continue
                else
                    for cntDatBlock=1:length(Ndatablock)
                        fpos=EnsStart{cntens}+...
                            DataOffset{cntens}(Ndatablock(cntDatBlock));
                        readNMEAVTG(fileid(cntens),fpos,cntens,cntDatBlock);
                    end
                end
            end
        end
        if any(AllHeaders(:,2)==WRII_DBT_ID)                               % Same as above but for WRII DBT data
            disp('WinRiver II DBT data found')
            dd=diff(find(AllHeaders(:,2)==WRII_DBT_ID &...
                AllHeaders(:,1)==8226));
            nblocks=max(diff(find(dd~=1)));
            initNMEADBT(nens,nblocks);
            for cntens=1:nens
                Ndatablock=find(DataHeader{cntens}(:,2)==WRII_DBT_ID &...
                    DataHeader{cntens}(:,1)==8226);
                if isempty(Ndatablock)
                    continue
                else
                    for cntDatBlock=1:length(Ndatablock)
                        fpos=EnsStart{cntens}+...
                            DataOffset{cntens}(Ndatablock(cntDatBlock));
                        readNMEADBT(fileid(cntens),fpos,cntens,cntDatBlock);
                    end
                end
            end
        end
    end
    if any(AllHeaders(:,1)==8448)                                           % Same as Velocity data but for WinRiver DBT data
        tmpdbt=cell(nens,1);
        for cntens=1:nens
            Ndatablock=find(DataHeader{cntens}(:,1)==8448);
            if isempty(Ndatablock) || ~helpers.load_nmea_package()
                continue
            else
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
                tmpdbt{cntens}=getNMEA(fileid(cntens),fpos);
            end
        end
        readNMEAint(tmpdbt,nmea.DBTMessage);
        clear tmpdbt;    
    end
    if any(AllHeaders(:,1)==8449)                                          % Sama as above but for WR GGA data
        tmpgga=cell(nens,1);
        for cntens=1:nens
            Ndatablock=find(DataHeader{cntens}(:,1)==8449);
            if isempty(Ndatablock) || ~helpers.load_nmea_package()
                continue
            else
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
                tmpgga{cntens}=getNMEA(fileid(cntens),fpos);
            end
        end
        readNMEAint(tmpgga,nmea.GGAMessage);
        clear tmpgga
    end
    if any(AllHeaders(:,1)==8450)                                          % Same as above but for WR VTG data
        tmpvtg=cell(nens,1);
        for cntens=1:nens
            Ndatablock=find(DataHeader{cntens}(:,1)==8450);
            if isempty(Ndatablock) || ~helpers.load_nmea_package()
                continue
            else
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
                tmpvtg{cntens}=getNMEA(fileid(cntens),fpos);
            end
        end
        readNMEAint(tmpvtg,nmea.VTGMessage);
        clear tmpvtg;
    end
    if any(AllHeaders(:,1)==8452)                                          % Same as above but for WR HDT data
        tmphdt=cell(nens,1);
        for cntens=1:nens
            Ndatablock=find(DataHeader{cntens}(:,1)==8452);
            if isempty(Ndatablock) || ~helpers.load_nmea_package()
                continue
            else
                fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
                tmphdt{cntens}=getNMEA(fileid(cntens),fpos);
            end
        end
        readNMEAint(tmphdt,nmea.HDTMessage);
        clear tmphdt;
    end
end


%% close all files
for cntfiles=1:nValidFiles
    fclose(ValidFilesId(cntfiles));
end


%% Output data
ADCPout=dataout;

clear global dataout WRII_GGA_ID WRII_HDT_ID WRII_VTG_ID WRII_DBT_ID

%% 
function valens=isvalens(fid,headpos)
% Check validity of ensemble

fmes=fseek(fid,headpos+2,-1);                                              %Go to byte after two header bytes
if fmes==-1                                                                %Check if move was possible
    valens=false;
    return
end
EnsBytes=fread(fid,1,'uint16=>double');                                    %read number of bytes in ensemble
if isempty(EnsBytes)
    valens=false;
    return
end
fmes=fseek(fid,headpos+EnsBytes+2,-1);                                     %Go to the end of the ensemble to check if file doesn't end before
if fmes==-1
    valens=false;
    return
end
fseek(fid,headpos,-1);                                                     %move to beginning of ensemble 
sumbytes=dec2bin(sum(fread(fid,EnsBytes,'uint8')));                        %calculate checksum
checksum=fread(fid,1,'*uint16');                                           %read checksum
if checksum==bin2dec(sumbytes(max((length(sumbytes)-15),1):end))                    %Check checksum
    valens=true;
else
    valens=false;
end

function headpos=searchHead(fid,fpos)
% search in the file for a valid header
headpos=-1;
fmes=fseek(fid,fpos,-1);
if fmes==-1
    return
end
while (~feof(fid))&& headpos==-1                                           %Loop until a end of file or until a header is found
    head=fread(fid,1,'*uint8');                                            %Read 1st byte
    if isempty(head) || head~=127                                          %Check for valid header value
        continue                                                           %If bad continue searching
    end
    head=fread(fid,1,'*uint8');                                            %Read 2nd byte
    if isempty(head)|| head~=127                                           %Check for valid header value
        continue                                                           %If bad continue searching
    end
    headpos=ftell(fid)-2;
end

%%          READING FUNCTIONS

% Function to check ensemble validity and read pos of ensemble and its data
function [NDataTypes,DataOffset,DataHeader]=readhead(fid,headpos)
fseek(fid,headpos+5,-1);                                                   %Move to the 'number of data-types'-field
NDataTypes=fread(fid,1,'uint8=>double');                                   %Read the number of data-types
DataOffset=fread(fid,NDataTypes,'uint16');                                 %Read offsets to data-types
DataHeader=ones(NDataTypes,2,'uint16')*65535;                              %Initialize header vector
for cntDataType=1:NDataTypes
    fseek(fid,headpos+DataOffset(cntDataType),-1);                         %Move to data header in ensemble
    DataHeader(cntDataType,:)=fread(fid,2,'*uint16')';                     %read header of data type (And next byte for navigation data)
end

% Function to read the fixed leader
function dataout = readFL(fid, fpos, cntens, dataout)
fseek(fid,fpos,-1);
LeadID=fread(fid,2,'*uint8');
if any(LeadID~=[0;0])
    error('readADCP:readFL:wrongID','Fixed Leader ID seems to be wrong')
end
dataout.firmver(cntens)=fread(fid,1,'*uint8');                                 %Firmware version
dataout.firmrev(cntens)=fread(fid,1,'*uint8');                                 %Firmware revision
dataout.sysconf(:,cntens)=num2str(fread(fid,16,'*ubit1'));                      %System configuration (see manual for explanation)
dataout.sysconfstr{cntens}=sysinterp(dataout.sysconf(:,cntens));                     %Interprete system configuration
dataout.SymData(cntens)=fread(fid,1,'*uint8');                                 %Flag for Real or Symulated data (0 for real data)
dataout.LagLength(cntens)=fread(fid,1,'*uint8');                               %Time period between sound pulses
dataout.usedbeams(cntens)=fread(fid,1,'*uint8');                               %Number of used beams
dataout.nbins(cntens)=fread(fid,1,'*uint8');                                   %number of bins (1-128)
dataout.pingperens(cntens)=fread(fid,1,'*uint16');                             %pings per ensemble (0-16384)
dataout.binsize(cntens)=fread(fid,1,'*uint16');                                %bin size in cm (1-6400)
dataout.blnk(cntens)=fread(fid,1,'*uint16');                                   %blanking in cm (0-9999)
fseek(fid,1,0);                                                            %skip data processing mode (always one)
dataout.minthrsh(cntens)=fread(fid,1,'*uint8');                                %Minimum threshold correlation in counts (1-256)
dataout.ncodrep(cntens)=fread(fid,1,'*uint8');                                 %Code repetitions in transmit pulse in counts (1-256)
dataout.minpercgood(cntens)=fread(fid,1,'*uint8');                             %Minimum percentage of good pings in one ensemble (1-100)
dataout.maxerrvel(cntens)=fread(fid,1,'*uint16');                              %Maximum value of error velocity in mm/s (1-5000 mm/s)
dataout.Tbetweenpng(cntens)=fread(fid,1,'uint8=>uint16')*6000+...
    fread(fid,1,'uint8=>uint16')*100+fread(fid,1,'uint8=>uint16');         %Time between two pings in cs (reading min, secs and cs)
dataout.corinfo(:,cntens)=num2str(fread(fid,8,'*ubit1'));                       %Coordinate information (see manual for explanation)
dataout.corstr{cntens}=corinterpr(dataout.corinfo(:,cntens));                        %Interprete coordinate information
dataout.headalign(cntens)=fread(fid,1,'*int16');                               %Head physical alignment correction in 0.01 degrees (-179.99 to 180.00)
dataout.headbias(cntens)=fread(fid,1,'*int16');                                %Head magnetic bias correction in 0.01 degrees (-179.99 to 180.00)
dataout.sensource(:,cntens)=num2str(fread(fid,8,'*ubit1'));                     %Sensor source information (see manual for explanation)
dataout.senavail(:,cntens)=num2str(fread(fid,8,'*ubit1'));                      %Sensor availability info  (see manual for explanation)
dataout.distmidbin1(cntens)=fread(fid,1,'*uint16');                            %Distance to middle of first bin in cm (0-65535)
dataout.lngthtranspulse(cntens)=fread(fid,1,'*uint16');                        %Length of the transmitted pulse in cm (0-65535)
dataout.watrefbins(:,cntens)=fread(fid,2,'*uint8');                              %Vector with begin and end bin for averaging to determine Water layer reference (1-128)
dataout.mintarget(cntens)=fread(fid,1,'*uint8');                               %Minimum for false target rejection in counts(0-255)
dataout.lowlattrig(cntens)=fread(fid,1,'*uint8');                              %Skip CX-command setting
dataout.distpulse(cntens)=fread(fid,1,'*uint16');                              %Distance between pulse repetitions in cm (0-65535) dependent on WM command
dataout.cpuserial(:,cntens)=fread(fid,8,'*uint8');                               %CPU serial number
dataout.bandwidth(cntens)=fread(fid,1,'*uint16');                              %Bandwidth (WB command)
dataout.syspower(cntens)=fread(fid,1,'*uint8');                                %System power in counts (CQ-command, only affects 75 and 150 KHz systems)
dataout.basefreqid(cntens)=fread(fid,1,'*uint8');                              %Base frequency index (only for Navigators)
dataout.serial(:,cntens)=fread(fid,4,'*uint8');                                  %ADCP serial number (REMUS only)
dataout.HADCPbeamangle(cntens)=fread(fid,1,'*uint8');                          %Beam angle, only for HADCP's

%Read velocity
function V=readVEL(fid,fpos,nbins)
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[0;1])
    error('readADCP:readVEL:wrongID','Velocity data ID seems to be wrong')
end
V=reshape(fread(fid,4*double(nbins),'*int16'),4,[])';

%Read echo intensity
function H=readECHO(fid,fpos,nbins)
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[0;3])
    error('readADCP:readECHO:wrongID','Echo intesity data ID seems to be wrong')
end
H=reshape(fread(fid,4*double(nbins),'*uint8'),4,[])';    

%Read correlation magnitude
function C=readCORR(fid,fpos,nbins)
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[0;2])
    error('readADCP:readCORR:wrongID','Correlation data ID seems to be wrong')
end
C=reshape(fread(fid,4*double(nbins),'*uint8'),4,[])';    

%Read percentage good
function P=readPERC(fid,fpos,nbins)
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[0;4])
    error('readADCP:readPERC:wrongID','Percentage good data ID seems to be wrong')
end
P=reshape(fread(fid,4*double(nbins),'*uint8'),4,[])';    

%Read ambient sound
function AS=readAS(fid,fpos)
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[2;12])% What should be --> just check by reading..... 
    error('readADCP:readAS:wrongID','Ambient sound data ID seems to be wrong')
end
AS=reshape(fread(fid,4*double(1),'*uint8'),4,[])';


% Read bottom track data
function readBT(fid,fpos,cntens)
global dataout
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[0;6])
    error('readADCP:readBT:wrongID','Bottom-Track data ID seems to be wrong')
end
bt1=fread(fid,2,'*uint16');
dataout.btpingperens(cntens)=bt1(1);                  %Bottom tracking pings per ensemble (0-999)
dataout.reacqdelay(cntens)=bt1(2);                    %Delay in number of ensembles before reacquiring (0-999)
bt2=fread(fid,4,'*uint8');
dataout.mincormag(cntens)=bt2(1);                      %Minimum for correlation magnitude in counts (0-255)
dataout.minevampl(cntens)=bt2(2);                      %Minimum evaluation amplitude in counts (1-255)
dataout.btminpergood(cntens)=bt2(3);                   %Minimum percentage of good bt pings 
dataout.btmode(cntens)=bt2(4);                         %BT mode
dataout.btmaxerrv(cntens)=fread(fid,1,'*uint16');                     %Maximum bt error velocity in mm/s (0-5000)
fseek(fid,4,0);                                                 %Reserved data
bt3=fread(fid,4,'uint16=>uint32');
dataout.btrange(cntens,1:4)=bt3(1:4);                      %bt range of beam 1,2,3,4 in cm (0-65535)
bt4=fread(fid,4,'*int16');
dataout.btvel(cntens,1:4)=bt4(1:4);                         %bt velocty of beam 1,2,3,4 in mm/s (-32768 to 32768)
bt5=fread(fid,12,'*uint8');
dataout.btcor(cntens,1:4)=bt5(1:4);                         %bt correlation magnitude beam 1,2,3,4 in counts (0-255)
dataout.btevampl(cntens,1:4)=bt5(5:8);                      %bt evaluation amplitude for strength of bottom echo beam 1,2,3,4 in counts (0-255)
dataout.btpercgood(cntens,1:4)=bt5(9:12);                    %bt percentage of good pings in beam 1,2,3,4 (0-100)
bt6=fread(fid,3,'*uint16');
dataout.minlyrsize(cntens)=bt6(1);                    %Minimum size of ref layer in dm (0-9999)
dataout.nearbnd(cntens)=bt6(2);                       %Near boundary of ref layer in dm (0-9999)
dataout.farbnd(cntens)=bt6(3);                        %Far boundary of ref layer in dm (0-9999) 
bt7=fread(fid,4,'*int16');    
dataout.reflyrvel(cntens,1:4)=bt7(1:4);                     %Reference layer velocity 1,2,3,4 in mm/s (-32768 to 32768)
bt8=fread(fid,12,'*uint8');
dataout.reflyrcor(cntens,1:4)=bt8(1:4);                     %Reference layer correlation magnitude beam 1,2,3,4 in counts (0-255)
dataout.reflyrint(cntens,1:4)=bt8(5:8);                     %Reference layer echo intensity beam 1,2,3,4 in counts (0-255)
dataout.reflyrpergood(cntens,1:4)=bt8(9:12);                 %Reference layer percentage good pings in beam 1,2,3,4 (0-100)
dataout.maxdepth(cntens)=fread(fid,1,'*uint16');                      %bt maximum depth in dm (80-9999)
bt9=fread(fid,5,'*uint8');
dataout.rssiamp(cntens,1:4)=bt9(1:4);                       %Received signal strength indicator in beam 1,2,3,4 in counts (0-255), 1 count ca. 0.45 dB
dataout.gain(cntens)=bt9(5);                           %Gain level for shallow water
bt10=fread(fid,4,'uint8=>uint32');
dataout.btrange(cntens,1:4)=reshape(bt10(1:4)*65636,[1,4])+dataout.btrange(cntens,1:4);     %Most significant byte of bt range in beam 1,2,3,4 in cm (65636-16777215)
if ~dataout.firmver(cntens) == rdi.ADCP_Type.STREAMPRO_31
    return % stop here if data is not from streampro
end
fseek(fid,4,0);
bt11=fread(fid,4,'*uint8');
dataout.sp_btrange_fract(cntens,1:4)=reshape(bt11,[1,4]);

% Read variable leader
function readVL(fid,fpos,cntens)
global dataout
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[128;0])
    error('readADCP:readVL:wrongID','Variable Leader ID seems to be wrong')
end
dataout.ensnum(cntens)=fread(fid,1,'uint16=>uint32');                                  %Ensemble number without rollover correction
% read first real time clock field without century information
vl7 = fread(fid,7,'uint8');
dataout.timeV1C(cntens,:)=[vl7(1) vl7(2) vl7(3) vl7(4) vl7(5) ...      
							    vl7(6)+vl7(7)/100];
dataout.ensnum(cntens)=fread(fid,1,'uint8=>uint32')*65535+dataout.ensnum(cntens);    %Ensemble number incuding rollover information
dataout.BITcheck(cntens,1:8)=num2str(fread(fid,8,'ubit1'))';                    %Result of Built in test (see manual for explanation)
fseek(fid,1,0);                                                            %Byte for BIT check which is reserved for future use
vl1=fread(fid,3,'*uint16');
dataout.speedsound(cntens)=vl1(1);                                              %Speed of sound in m/s (1400-1600)
dataout.depthtransd(cntens)=vl1(2);                                             %Depth of transducer in decimeters (1-9999)
dataout.heading(cntens)=vl1(3);                                                 %Heading in 0.01 degrees (000.00-359.99)
vl2=fread(fid,2,'*int16');
dataout.pitch(cntens)=vl2(1);                                                   %Pitch in 0.01 degrees (-20.00 to 20.00)
dataout.roll(cntens)=vl2(2);                                                    %Roll in 0.01 degrees (-20.00 to 20.00)
dataout.salinity(cntens)=fread(fid,1,'*uint16');                                %Salinity in parts per thousands (0-40)
dataout.temperature(cntens)=fread(fid,1,'*int16');                              %Temeperature in 0.01 degrees celcius (-5.00 to 40.00)
vl3=fread(fid,14,'*uint8');
dataout.prepingT(cntens)=uint16(vl3(1))*600+uint16(vl3(2))*100+uint16(vl3(3));  %Pre ping time in cs
dataout.headstd(cntens)=vl3(4);                                                 %Standard deviation in heading in degrees (0-180)
dataout.pitchstd(cntens)=vl3(5);                                                %Standard deviation in pitch in 0.1 degrees (0.0-20.0)
dataout.rollstd(cntens)=vl3(6);                                                 %Standard deviation in roll in 0.1 degrees (0.0-20.0)
dataout.ADC(cntens,1:8)=vl3(7:14);                                              %ADC channels, each columns is one channel
vl4=fread(fid,32,'*ubit1');
dataout.errorstat1(cntens,1:8)=vl4(1:8);                                        %Error status check (for explanation see manual)
dataout.errorstat2(cntens,1:8)=vl4(9:16);                                       %Error status check (for explanation see manual)
dataout.errorstat3(cntens,1:8)=vl4(17:24);                                      %Error status check (for explanation see manual)
dataout.errorstat4(cntens,1:8)=vl4(25:32);                                      %Error status check (for explanation see manual)
fseek(fid,2,0);                                                                 %Reserved for RDI
vl5=fread(fid,2,'*uint32');
if dataout.firmver(cntens) == rdi.ADCP_Type.STREAMPRO_31
    return % for strempro, data ends here
end
dataout.pressure(cntens)=vl5(1);                                                %Pressure in decapascal (0-4'294'967'295)
dataout.pressurevar(cntens)=vl5(2);                                             %Variance in pressure in decapascal (0-4'294'967'295)
fseek(fid,1,0);                                                            %Spare byte
vl6=fread(fid,8,'uint8');
dataout.timeV(cntens,:)=[vl6(1)*100+vl6(2) vl6(3) vl6(4) vl6(5) vl6(6) ...
    vl6(7)+vl6(8)/100];                                                    %Gives date in Matlab vector

function readSP(fid,fpos,cntens)
global dataout
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[0;80])
    error('readADCP:readSP:wrongID','StreamPro Leader ID seems to be wrong')
end
dataout.sp_longlag(cntens)=fread(fid,1,'*uint16');
dataout.sp_shortlag(cntens)=fread(fid,1,'*uint16');
dataout.sp_pgood(cntens)=fread(fid,1,'*uint16');
dataout.sp_subpings(cntens)=fread(fid,1,'*uint16');
dataout.sp_lastbin(cntens)=fread(fid,1,'*uint16');
dataout.sp_corr_thresh(cntens)=fread(fid,1,'*uint8');
dataout.sp_mid_bin1(cntens)=fread(fid,1,'*uint16');
dataout.sp_bin_size(cntens)=fread(fid,1,'*uint16');
dataout.sp_bin_space(cntens)=fread(fid,1,'*uint16');
dataout.sp_transmit_length(cntens)=fread(fid,1,'*uint16');

function readTM(fid,fpos,cntens)
global dataout
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[0;50])
    error('readADCP:readTM:wrongID','Transformation matrix ID seems to be wrong')
end
dataout.transformation_matrix(:,:,cntens)=reshape(fread(fid,16,'*int16'),4,4);

function readVBeamRange(fid,fpos,cntens)
    global dataout
    fseek(fid,fpos+2,-1);
    dataout.vbeam_eval_amp(cntens) = fread(fid,1,'*uint8');
    dataout.vbeam_rssi_amp(cntens) = fread(fid,1,'*uint8');
    dataout.vbeam_range(cntens) = fread(fid,1,'*uint32');
    dataout.vbeam_status(cntens) = fread(fid,1,'*uint8');

function readVBeamLeader(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
dataout.vbeam_ncells(cntens) = fread(fid, 1, '*uint16');
dataout.vbeam_ping_per_ens(cntens) = fread(fid, 1, '*uint16');
dataout.vbeam_cellsize = fread(fid, 1, '*uint16');
dataout.vbeam_distmidbin1(cntens) = fread(fid, 1, '*uint16');
fseek(fid,2,0);
dataout.vbeam_xmit_length(cntens) = fread(fid, 1, '*uint16');
dataout.vbeam_lag_length(cntens) = fread(fid,1,'*uint16');
dataout.vbeam_ncode_xmit(cntens) = fread(fid, 1, '*uint16');

function readVBeamVelocity(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
ncells = dataout.vbeam_ncells(cntens);
dataout.vbeam_velocity(1:ncells,cntens) = ...
    fread(fid,ncells,"*int16");

function readVBeamCorr(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
ncells = dataout.vbeam_ncells(cntens);
dataout.vbeam_corr(1:ncells,cntens) = ...
    fread(fid,ncells,"*uint8");

function readVBeamEcho(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
ncells = dataout.vbeam_ncells(cntens);
dataout.vbeam_echo(1:ncells,cntens) = ...
    fread(fid,ncells,"*uint8");

function readVBeamPerc(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
ncells = dataout.vbeam_ncells(cntens);
dataout.vbeam_perc(1:ncells,cntens) = ...
    fread(fid,ncells,"*uint8");

function readVBeamStat(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
ncells = dataout.vbeam_ncells(cntens);
dataout.vbeam_stat(1:ncells,cntens) = ...
    fread(fid,ncells,"*uint8");

function readSLayerVelLeader(fid,fpos,cntens)
    global dataout
    fseek(fid,fpos+2,-1);
    dataout.sl_ncells(cntens) = fread(fid,1,'*uint8');
    dataout.sl_cellsize(cntens) = fread(fid,1,'*uint16');
    dataout.sl_distmidbin1(cntens) = fread(fid,1,'*uint16');

function readSLayerVelocity(fid,fpos,cntens)
    global dataout
    fseek(fid,fpos+2,-1);
    ncells = dataout.sl_ncells(cntens);
    nvels = 4 * dataout.sl_ncells(cntens);
    dataout.sl_velocity(1:ncells,cntens,:) = permute(reshape(...
        fread(fid,nvels,"*int16"),4,ncells),[2 3 1]);

function readSLayerEcho(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
ncells = dataout.sl_ncells(cntens);
nvels = 4 * dataout.sl_ncells(cntens);
dataout.sl_echo(1:ncells,cntens,:) = permute(reshape(...
    fread(fid,nvels,"*uint8"),4,ncells),[2 3 1]);

function readSLayerCorr(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
ncells = dataout.sl_ncells(cntens);
nvels = 4 * dataout.sl_ncells(cntens);
dataout.sl_corr(1:ncells,cntens,:) = permute(reshape(...
    fread(fid,nvels,"*uint8"),4,ncells),[2 3 1]);

function readSLayerPerc(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
ncells = dataout.sl_ncells(cntens);
nvels = 4 * dataout.sl_ncells(cntens);
dataout.sl_perc(1:ncells,cntens,:) = permute(reshape(...
    fread(fid,nvels,"*uint8"),4,ncells),[2 3 1]);

function readSLayerStat(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
ncells = dataout.sl_ncells(cntens);
nvels = 4 * dataout.sl_ncells(cntens);
dataout.sl_stat(1:ncells,cntens,:) = permute(reshape(...
    fread(fid,nvels,"*uint8"),4,ncells),[2 3 1]);

function readAutoMode3(fid,fpos,cntens)
global dataout
fseek(fid,fpos+2,-1);
dataout.am3_nbeams(cntens) = fread(fid,1,'uint8');
for cb  = 1:dataout.am3_nbeams(cntens)
    dataout.am3_setup(1,cntens,cb) = fread(fid,1,'*uint8');
    dataout.am3_depth(1,cntens,cb) = fread(fid,1,'*uint16');
    dataout.am3_ping_count(1,cntens,cb) = fread(fid,1,'*uint8');
    dataout.am3_ping_type(1,cntens,cb) = fread(fid,1,'*uint8');
    dataout.am3_ncells(1,cntens,cb) = fread(fid,1,'*uint16');
    dataout.am3_cellsize(1,cntens,cb) = fread(fid,1,'*uint16');
    dataout.am3_distmidbin1(1,cntens,cb) = fread(fid,1,'*uint16');
    dataout.am3_ncodereps(1,cntens,cb) = fread(fid,1,'*uint8');
    dataout.am3_xmit_length(1,cntens,cb) = fread(fid,1,'*uint16');
    dataout.am3_lag_length(1,cntens,cb) = fread(fid,1,'*uint16');
    dataout.am3_xmit_bandwidth(1,cntens,cb) = fread(fid,1,'*uint8');
    dataout.am3_recv_bandwidth(1,cntens,cb) = fread(fid,1,'*uint8');
    dataout.am3_min_ping_interval(1,cntens,cb) = fread(fid,1,'*uint16');
end

% Read WinRiverII General NMEA GGA data
function readNMEAGGA(fid,fpos,cntens,cntblock)
global dataout WRII_GGA_ID
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint16');
if any(dataID~=[8226;WRII_GGA_ID])
    error('readADCP:readNMEAGGA:wrongID','WinRiver II general NMEA GGA data ID seems to be wrong')
end
fseek(fid,2,0);
dataout.NMEAGGA.deltaT(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEAGGA.msgHeader{cntens}(cntblock,:)=fread(fid,7,'uchar=>char')';
dataout.NMEAGGA.UTC{cntens}(cntblock,:)=fread(fid,10,'uchar=>char')';
dataout.NMEAGGA.Lat(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEAGGA.SN(cntens,cntblock)=fread(fid,1,'uchar=>char');
dataout.NMEAGGA.Long(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEAGGA.EW(cntens,cntblock)=fread(fid,1,'uchar=>char');
dataout.NMEAGGA.Qual(cntens,cntblock)=fread(fid,1,'*uint8');
dataout.NMEAGGA.NSat(cntens,cntblock)=fread(fid,1,'*uint8');
dataout.NMEAGGA.HDOP(cntens,cntblock)=fread(fid,1,'*float32');
dataout.NMEAGGA.Alt(cntens,cntblock)=fread(fid,1,'*float32');
dataout.NMEAGGA.AltUnit(cntens,cntblock)=fread(fid,1,'uchar=>char');
dataout.NMEAGGA.Geoid(cntens,cntblock)=fread(fid,1,'*float32');
dataout.NMEAGGA.GeoidUnit(cntens,cntblock)=fread(fid,1,'uchar=>char');
dataout.NMEAGGA.AgeDGPS(cntens,cntblock)=fread(fid,1,'*float32');
dataout.NMEAGGA.RefStatID(cntens,cntblock)=fread(fid,1,'*uint16');

% Read WinRiverII General NMEA HDT data
function readNMEAHDT(fid,fpos,cntens,cntblock)
global dataout WRII_HDT_ID
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint16');
if any(dataID~=[8226;WRII_HDT_ID])
    error('readADCP:readNMEAHDT:wrongID','WinRiver II general NMEA HDT data ID seems to be wrong')
end
fseek(fid,2,0);
dataout.NMEAHDT.deltaT(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEAHDT.msgHeader{cntens}(cntblock,:)=fread(fid,7,'uchar=>char')';
dataout.NMEAHDT.heading(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEAHDT.trueInd(cntens,cntblock)=fread(fid,1,'uchar=>char');

% Read WinRiverII General NMEA VTG data
function readNMEAVTG(fid,fpos,cntens,cntblock)
global dataout WRII_VTG_ID
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint16');
if any(dataID~=[8226;WRII_VTG_ID])
    error('readADCP:readNMEAVTG:wrongID','WinRiver II general NMEA VTG data ID seems to be wrong')
end
fseek(fid,2,0);
dataout.NMEAVTG.deltaT(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEAVTG.msgHeader{cntens}(cntblock,:)=fread(fid,7,'uchar=>char')';
dataout.NMEAVTG.COGTrue(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEAVTG.TrueIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
dataout.NMEAVTG.COGMagn(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEAVTG.MagnIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
dataout.NMEAVTG.SpeedOverGroundKts(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEAVTG.KtsIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
dataout.NMEAVTG.SpeedOverGroundKmh(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEAVTG.KmhIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
dataout.NMEAVTG.ModeIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');

% Read WinRiverII General NMEA DBT data
function readNMEADBT(fid,fpos,cntens,cntblock)
global dataout WRII_DBT_ID
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint16');
if any(dataID~=[8226;WRII_DBT_ID])
    error('readADCP:readNMEADBT:wrongID','WinRiver II general NMEA DBT data ID seems to be wrong')
end
fseek(fid,2,0);
dataout.NMEADBT.deltaT(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEADBT.msgHeader{cntens}(cntblock,:)=fread(fid,7,'uchar=>char')';
dataout.NMEADBT.WaterDepthFt(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEADBT.FeetIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
dataout.NMEADBT.WaterDepthm(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEADBT.MeterIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
dataout.NMEADBT.WaterDepthF(cntens,cntblock)=fread(fid,1,'*float64');
dataout.NMEADBT.FathomIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');

% Read WinRiver GGA data
function ggastr=getGGA(fid,fpos)
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[1;33])
    error('readADCP:getGGA:wrongID','WinRiver GGA data ID seems to be wrong')
end
ggastr=fgetl(fid);

function dat=expand_nmea(has_data, dat)
data_fields = fieldnames(dat);
for cm = 1:numel(data_fields)
    cur_dat = dat.(data_fields{cm});
    if isfloat(cur_dat)
        dat.(data_fields{cm}) = nan(size(cur_dat));
    end
    dat.(data_fields{cm})(has_data,:) = cur_dat;
end

function readNMEAint(ggastr,msg)
global dataout
hasdata=~cellfun(@isempty,regexp(ggastr,msg.message_pattern));
if ~any(hasdata)
    return
end
datgga=msg.parse([ggastr{hasdata}]);
if ~isfield(datgga, msg.name)
    return
end
datgga.(msg.name) = expand_nmea(hasdata,datgga.(msg.name));
dataout.(msg.name) = datgga.(msg.name);

%Read WinRiver HDT data
function nmeastr=getNMEA(fid,fpos)
fseek(fid,fpos+2,-1);
nmeastr=fgetl(fid);


%Read WinRiver DBT data
function dbtstr=getDBT(fid,fpos)
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[0;33])
    error('readADCP:getDBT:wrongID','WinRiver DBT data ID seems to be wrong')
end
dbtstr=fgetl(fid);

function readDBTint(dbtstr,nens)
global dataout
rdi.defineNMEA;
hasdata=~cellfun(@isempty,regexp(dbtstr,patterns.dbt));
if ~any(hasdata)
    return
end
disp('WinRiver DBT data found')
initDBT(nens);
datdbt=rdi.readDBT(dbtstr(hasdata));
dataout.DBT.depthf(hasdata)=datdbt.depthf;
dataout.DBT.dephtM(hasdata)=datdbt.depthM;
dataout.DBT.depthF(hasdata)=datdbt.depthF;


%Read WinRiver VTG data
function vtgstr=getVTG(fid,fpos)
fseek(fid,fpos,-1);
dataID=fread(fid,2,'*uint8');
if any(dataID~=[2;33])
    error('readADCP:getVTG:wrongID','WinRiver VTG data ID seems to be wrong')
end
vtgstr=fgetl(fid);

function readVTGint(vtgstr,nens)
global dataout
rdi.defineNMEA
hasdata=~cellfun(@isempty,regexp(vtgstr,patterns.vtg));
if ~any(hasdata)
    return
end
disp('WinRiver VTG data found')
initVTG(nens);
datvtg=rdi.readVTG(vtgstr(hasdata));
dataout.VTG.TrackDegTrue(hasdata)=datvtg.TrackDegTrue;
dataout.VTG.TrackDegMagn(hasdata)=datvtg.TrackDegMagn;
dataout.VTG.SpeedKnots(hasdata)=datvtg.SpeedKnots;
dataout.VTG.SpeedKmH(hasdata)=datvtg.SpeedKmH;
dataout.VTG.mode(hasdata)=datvtg.mode;


%%          PREALLOCATION FUNCTIONS

function dataout = initFL(nens, dataout)
dataout.firmver=zeros(1,nens,'uint8');                                 %Firmware version
dataout.firmrev=zeros(1,nens,'uint8');                                 %Firmware revision
dataout.sysconf=repmat('0',[16, nens]);                      %System configuration (see manual for explanation)
dataout.sysconfstr=cell(1,nens);                     %Interprete system configuration
dataout.SymData=zeros(1,nens,'uint8');                                 %Flag for Real or Symulated data (0 for real data)
dataout.LagLength=zeros(1,nens,'uint8');                               %Time period between sound pulses
dataout.usedbeams=zeros(1,nens,'uint8');                               %Number of used beams
dataout.nbins=zeros(1,nens,'uint8');                                   %number of bins (1-128)
dataout.pingperens=zeros(1,nens,'uint16');                             %pings per ensemble (0-16384)
dataout.binsize=zeros(1,nens,'uint16');                                %bin size in cm (1-6400)
dataout.blnk=zeros(1,nens,'uint16');                                   %blanking in cm (0-9999)
dataout.minthrsh=zeros(1,nens,'uint8');                                %Minimum threshold correlation in counts (1-256)
dataout.ncodrep=zeros(1,nens,'uint8');                                 %Code repetitions in transmit pulse in counts (1-256)
dataout.minpercgood=zeros(1,nens,'uint8');                             %Minimum percentage of good pings in one ensemble (1-100)
dataout.maxerrvel=zeros(1,nens,'uint16');                              %Maximum value of error velocity in mm/s (1-5000 mm/s)
dataout.Tbetweenpng=zeros(1,nens,'uint16');
dataout.corinfo=repmat('0',[8, nens]);                       %Coordinate information (see manual for explanation)
dataout.corstr=cell(1,nens);
dataout.headalign=zeros(1,nens,'int16');                               %Head physical alignment correction in 0.01 degrees (-179.99 to 180.00)
dataout.headbias=zeros(1,nens,'int16');                                %Head magnetic bias correction in 0.01 degrees (-179.99 to 180.00)
dataout.sensource=repmat('0',[8, nens]);                     %Sensor source information (see manual for explanation)
dataout.senavail=repmat('0',[8, nens]);                      %Sensor availability info  (see manual for explanation)
dataout.distmidbin1=zeros(1,nens,'uint16');                            %Distance to middle of first bin in cm (0-65535)
dataout.lngthtranspulse=zeros(1,nens,'uint16');                        %Length of the transmitted pulse in cm (0-65535)
dataout.watrefbins=zeros(2,nens,'uint8');                              %Vector with begin and end bin for averaging to determine Water layer reference (1-128)
dataout.mintarget=zeros(1,nens,'uint8');                               %Minimum for false target rejection in counts(0-255)
dataout.lowlattrig=zeros(1,nens,'uint8');                              %Skip CX-command setting
dataout.distpulse=zeros(1,nens,'uint16');                              %Distance between pulse repetitions in cm (0-65535) dependent on WM command
dataout.cpuserial=zeros(8,nens,'uint8');                               %CPU serial number
dataout.bandwidth=zeros(1,nens,'uint16');                              %Bandwidth (WB command)
dataout.syspower=zeros(1,nens,'uint8');                                %System power in counts (CQ-command, only affects 75 and 150 KHz systems)
dataout.basefreqid=zeros(1,nens,'uint8');                              %Base frequency index (only for Navigators)
dataout.serial=zeros(4,nens,'uint8');                                  %ADCP serial number (REMUS only)
dataout.HADCPbeamangle=zeros(1,nens,'uint8');                          %Beam angle, only for HADCP's



%Initialize Variable leader
function initVL(nens)
global dataout
dataout.ensnum= zeros(1,nens, 'uint32');
dataout.BITcheck= char(zeros(nens,8));
dataout.speedsound= zeros(1,nens, 'uint16');
dataout.depthtransd= zeros(1,nens, 'uint16');
dataout.heading= zeros(1,nens, 'uint16');
dataout.pitch= zeros(1,nens, 'int16');
dataout.roll= zeros(1,nens, 'int16');
dataout.salinity= zeros(1,nens, 'uint16');
dataout.temperature= zeros(1,nens, 'int16');
dataout.prepingT= zeros(1,nens, 'uint16');
dataout.headstd= zeros(1,nens, 'uint8');
dataout.pitchstd= zeros(1,nens, 'uint8');
dataout.rollstd= zeros(1,nens, 'uint8');
dataout.ADC= zeros(nens,8, 'uint8');
dataout.errorstat1= zeros(nens,8, 'uint8');
dataout.errorstat2= zeros(nens,8, 'uint8');
dataout.errorstat3= zeros(nens,8, 'uint8');
dataout.errorstat4= zeros(nens,8, 'uint8');
dataout.pressure= zeros(1,nens, 'uint32');
dataout.pressurevar= zeros(1,nens, 'uint32');
dataout.timeV= zeros(nens,6, 'double');

%Initialize StreamPro leader data
function initSP(nens)
global dataout
dataout.sp_longlag=zeros(1,nens, 'uint16');
dataout.sp_shortlag=zeros(1,nens, 'uint16');
dataout.sp_pgood=zeros(1,nens, 'uint16');
dataout.sp_subpings=zeros(1,nens, 'uint16');
dataout.sp_lastbin=zeros(1,nens, 'uint16');
dataout.sp_corr_thresh=zeros(1,nens, 'uint8');
dataout.sp_mid_bin1=zeros(1,nens, 'uint16');
dataout.sp_bin_size=zeros(1,nens, 'uint16');
dataout.sp_bin_space=zeros(1,nens, 'uint16');
dataout.sp_transmit_length=zeros(1,nens, 'uint16');

% Initialize transformation matrix data
function initTM(nens)
global dataout
dataout.transformation_matrix=zeros(4,4,nens,'int16');

function initVBeamRange(nens)
global dataout
dataout.vbeam_eval_amp = zeros(1,nens,'uint8');
dataout.vbeam_rssi_amp = zeros(1,nens,'uint8');
dataout.vbeam_range = zeros(1,nens,'uint32');
dataout.vbeam_status = zeros(1,nens,'uint8');

function initVBeamLeader(nens)
global dataout
dataout.vbeam_ncells = zeros(1,nens,'uint16');
dataout.vbeam_ping_per_ens = zeros(1,nens,'uint16');
dataout.vbeam_cell_size = zeros(1,nens,'uint16');
dataout.vbeam_distmidbin1 = zeros(1,nens,'uint16');
dataout.vbeam_xmit_length = zeros(1,nens,'uint16');
dataout.vbeam_lag_length = zeros(1,nens,'uint16');
dataout.vbeam_ncode_xmit = zeros(1,nens,'uint16');

function initVBeamVelocity(nens)
global dataout
ncells = max(dataout.vbeam_ncells);
dataout.vbeam_velocity = ones(ncells,nens,"int16")*intmin("int16");

function initVBeamPerc(nens)
global dataout
ncells = max(dataout.vbeam_ncells);
dataout.vbeam_perc = zeros(ncells,nens,"uint8");

function initVBeamEcho(nens)
global dataout
ncells = max(dataout.vbeam_ncells);
dataout.vbeam_echo = zeros(ncells,nens,"uint8");

function initVBeamCorr(nens)
global dataout
ncells = max(dataout.vbeam_ncells);
dataout.vbeam_corr = zeros(ncells,nens,"uint8");

function initVBeamStat(nens)
global dataout
ncells = max(dataout.vbeam_ncells);
dataout.vbeam_stat = zeros(ncells,nens,"uint8");

function initSLayerVelLeader(nens)
global dataout
dataout.sl_ncells = zeros(1,nens,'uint8');
dataout.sl_cellsize = zeros(1,nens,'uint16');
dataout.sl_distmidbin1 = zeros(1,nens,'uint16');

function initSLayerVelocity(nens)
global dataout
ncells = max(dataout.sl_ncells);
dataout.sl_velocity = ones(ncells,nens,4,"int16")*intmin("int16");

function initSLayerEcho(nens)
global dataout
ncells = max(dataout.sl_ncells);
dataout.sl_echo = zeros(ncells,nens,4,"uint8");

function initSLayerCorr(nens)
global dataout
ncells = max(dataout.sl_ncells);
dataout.sl_corr = zeros(ncells,nens,4,"uint8");

function initSLayerPerc(nens)
global dataout
ncells = max(dataout.sl_ncells);
dataout.sl_perc = zeros(ncells,nens,4,"uint8");

function initSLayerStat(nens)
global dataout
ncells = max(dataout.sl_ncells);
dataout.sl_stat = zeros(ncells,nens,4,"uint8");

function initAutoMode3(nens)
global dataout
dataout.am3_nbeams = zeros(1,nens,'uint8');
dataout.am3_setup = zeros(1,nens,4,'uint8');
dataout.am3_depth = zeros(1,nens,4,'uint16');
dataout.am3_ping_count = zeros(1,nens,4,'uint8');
dataout.am3_ping_type = zeros(1,nens,4,'uint8');
dataout.am3_ncells = zeros(1,nens,4,'uint16');
dataout.am3_cellsize = zeros(1,nens,4,'uint16');
dataout.am3_distmidbin1 = zeros(1,nens,4,'uint16');
dataout.am3_ncodereps = zeros(1,nens,4,'uint8');
dataout.am3_xmit_length = zeros(1,nens,4,'uint16');
dataout.am3_lag_length = zeros(1,nens,4,'uint16');
dataout.am3_xmit_bandwidth = zeros(1,nens,4,'uint8');
dataout.am3_recv_bandwidth = zeros(1,nens,4,'uint8');
dataout.am3_min_ping_interval = zeros(1,nens,4,'uint16');


%Initialize Bottom tracking data
function initBT(nens)
global dataout
dataout.btpingperens=zeros(1,nens,'uint16');
dataout.reacqdelay=zeros(1,nens,'uint16');
dataout.mincormag=zeros(1,nens,'uint8');
dataout.minevampl=zeros(1,nens,'uint8');
dataout.btminpergood=zeros(1,nens,'uint8');
dataout.btmode=zeros(1,nens,'uint8');
dataout.btmaxerrv=zeros(1,nens,'uint16');
dataout.btrange=zeros([nens,4],'uint32');
dataout.btvel=zeros([nens,4],'int16');
dataout.btcor=zeros([nens,4],'uint8');
dataout.btevampl=zeros([nens,4],'uint8');
dataout.btpercgood=zeros([nens,4],'uint8');
dataout.minlyrsize=zeros(1,nens,'uint16');
dataout.nearbnd=zeros(1,nens,'uint16');
dataout.farbnd=zeros(1,nens,'uint16');
dataout.reflyrvel=zeros([nens,4],'int16');
dataout.reflyrcor=zeros([nens,4],'uint8');
dataout.reflyrint=zeros([nens,4],'uint8');
dataout.reflyrpergood=zeros([nens,4],'uint8');
dataout.maxdepth=zeros(1,nens,'uint16');
dataout.rssiamp=zeros([nens,4],'uint8');
dataout.gain=zeros(1,nens,'uint8');
dataout.sp_btrange_fract=zeros([nens,4],'uint8');

%Initialize WinRiverII external GGA data
function initNMEAGGA(nens,nblocks)
global dataout
dataout.NMEAGGA.deltaT=ones(nens,nblocks,'double')*NaN;
dataout.NMEAGGA.msgHeader=cell(nens,1);      %char(ones(nens,7,'uint8'));
dataout.NMEAGGA.UTC=cell(nens,1);            %char(ones(nens,10,'uint8'));
dataout.NMEAGGA.Lat=ones(nens,nblocks,'double')*NaN;
dataout.NMEAGGA.SN=char(ones(nens,nblocks,'uint8'));
dataout.NMEAGGA.Long=ones(nens,nblocks,'double')*NaN;
dataout.NMEAGGA.EW=char(ones(nens,nblocks,'uint8'));
dataout.NMEAGGA.Qual=ones(nens,nblocks,'uint8');
dataout.NMEAGGA.NSat=ones(nens,nblocks,'uint8');
dataout.NMEAGGA.HDOP=ones(nens,nblocks,'single');
dataout.NMEAGGA.Alt=ones(nens,nblocks,'single')*NaN;
dataout.NMEAGGA.AltUnit=char(ones(nens,nblocks,'uint8'));
dataout.NMEAGGA.Geoid=ones(nens,nblocks,'single')*NaN;
dataout.NMEAGGA.GeoidUnit=char(ones(nens,nblocks,'uint8'));
dataout.NMEAGGA.AgeDGPS=ones(nens,nblocks,'single')*NaN;
dataout.NMEAGGA.RefStatID=ones(nens,nblocks,'uint16');

%Initialize WinRiverII external GGA data
function initNMEAHDT(nens,nblocks)
global dataout
dataout.NMEAHDT.deltaT=ones(nens,nblocks,'double')*NaN;
dataout.NMEAHDT.msgHeader=cell(nens,1);
dataout.NMEAHDT.heading=ones(nens,nblocks,'double')*NaN;
dataout.NMEAHDT.trueInd=char(ones(nens,nblocks,'uint8'));

%Initialize WinRiverII external GGA data
function initNMEAVTG(nens,nblocks)
global dataout
dataout.NMEAVTG.deltaT=ones(nens,nblocks,'double')*NaN;
dataout.NMEAVTG.msgHeader=cell(nens,1);
dataout.NMEAVTG.COGTrue=ones(nens,nblocks,'double')*NaN;
dataout.NMEAVTG.TrueIndicator=char(ones(nens,nblocks,'uint8'));
dataout.NMEAVTG.COGMagn=ones(nens,nblocks,'double')*NaN;
dataout.NMEAVTG.MagnIndicator=char(ones(nens,nblocks,'uint8'));
dataout.NMEAVTG.SpeedOverGroundKts=ones(nens,nblocks,'double')*NaN;
dataout.NMEAVTG.KtsIndicator=char(ones(nens,nblocks,'uint8'));
dataout.NMEAVTG.SpeedOverGroundKmh=ones(nens,nblocks,'double')*NaN;
dataout.NMEAVTG.KmhIndicator=char(ones(nens,nblocks,'uint8'));
dataout.NMEAVTG.ModeIndicator=char(ones(nens,nblocks,'uint8'));

%Initialize WinRiverII external GGA data
function initNMEADBT(nens,nblocks)
global dataout
dataout.NMEADBT.deltaT=ones(nens,nblocks,'double')*NaN;
dataout.NMEADBT.msgHeader=cell(nens,1);
dataout.NMEADBT.WaterDepthFt=ones(nens,nblocks,'double')*NaN;
dataout.NMEADBT.FeetIndicator=char(ones(nens,nblocks,'uint8'));
dataout.NMEADBT.WaterDepthm=ones(nens,nblocks,'double')*NaN;
dataout.NMEADBT.MeterIndicator=char(ones(nens,nblocks,'uint8'));
dataout.NMEADBT.WaterDepthF=ones(nens,nblocks,'double')*NaN;
dataout.NMEADBT.FathomIndicator=char(ones(nens,nblocks,'uint8'));

%Initialize WinRiver DBT external data
function initDBT(nens)
global dataout
dataout.DBT.depthf=nan(1,nens,'single');
dataout.DBT.dephtM=nan(1,nens,'single');
dataout.DBT.depthF=nan(1,nens,'single');

%Initialize WinRiver GGA external data
function initGGA(nens)
global dataout
dataout.GGA.UTCtime=nan(nens,3,'single');
dataout.GGA.lat=nan(1,nens,'double');
dataout.GGA.long=nan(1,nens,'double');
dataout.GGA.qualind=zeros(1,nens,'uint8');
dataout.GGA.numsat=zeros(1,nens,'uint8');
dataout.GGA.hdop=nan(1,nens,'single');
dataout.GGA.antalt=nan(1,nens,'single');
dataout.GGA.geosep=nan(1,nens,'single');
dataout.GGA.agediff=nan(1,nens,'single');
dataout.GGA.diffid=zeros(1,nens,'uint16');

%Initialize WinRiver VTG external data
function initVTG(nens)
global dataout
dataout.VTG.TrackDegTrue=nan(1,nens,'single');
dataout.VTG.TrackDegMagn=nan(1,nens,'single');
dataout.VTG.SpeedKnots=nan(1,nens,'single');
dataout.VTG.SpeedKmH=nan(1,nens,'single');
dataout.VTG.mode=zeros(1,nens,'uint8');

%Initialize WinRiver HDT external data
function initHDT(nens)
global dataout
dataout.HDT.heading=nan(1,nens,'single');


%%          INTERPRATATION FUNCTIONS

% Function to interpret system id
function sysstr=sysinterp(sysid)
    sysstrtemp='';
    switch sysid(1:3)'
        case '000'
            sysstrtemp=[sysstrtemp,'75 KHz System'];
        case '100'
            sysstrtemp=[sysstrtemp,'150 KHz System'];
        case '010'
            sysstrtemp=[sysstrtemp,'300 KHz System'];
        case '110'
            sysstrtemp=[sysstrtemp,'600 KHz System'];
        case '001'
            sysstrtemp=[sysstrtemp,'1200 KHz System'];
        case '101'
            sysstrtemp=[sysstrtemp,'2400 KHz System'];
    end
    if sysid(4)=='0'
        sysstrtemp=[sysstrtemp,', Concave beam pattern'];
    else
        sysstrtemp=[sysstrtemp,', Convex beam pattern'];
    end
    switch sysid(5:6)'
        case '00'
        sysstrtemp=[sysstrtemp,', Sensor configuration #1'];
        case '10'
        sysstrtemp=[sysstrtemp,', Sensor configuration #2'];
        case '11'
        sysstrtemp=[sysstrtemp,', Sensor configuration #3'];
    end
    if sysid(7)=='0'
        sysstrtemp=[sysstrtemp,', XDCR HD Not attached'];
    else
        sysstrtemp=[sysstrtemp,', XDCR HD Attached'];
    end
    if sysid(8)=='0'
        sysstrtemp=[sysstrtemp,', Up facing beam'];
    else
        sysstrtemp=[sysstrtemp,', Down facing beam'];
    end
    switch sysid(9:10)'
        case '00'
        sysstrtemp=[sysstrtemp,', Beam angle: 15E'];
        case '10'
        sysstrtemp=[sysstrtemp,', Beam angle: 20E'];
        case '11'
        sysstrtemp=[sysstrtemp,', Beam angle: 30E'];
        case '01'
        sysstrtemp=[sysstrtemp,', Unknown beam angle'];
    end
    switch sysid(13:16)'
        case '0010'
            sysstrtemp=[sysstrtemp,', 4 Beam Janus configuration'];
        case '1010'
            sysstrtemp=[sysstrtemp,', 5 Beam Janus configuration, 3 Demods'];
        case '1111'
            sysstrtemp=[sysstrtemp,', 5 Beam Janus configuration, 2 Demods'];
    end
sysstr=sysstrtemp;

% Function to interprete coordinate information (bits are reversed wrt manual of ADCP)
function corstr=corinterpr(corid)
            corstrtemp='';
            switch corid(4:5)'
                case '00'
                    corstrtemp=[corstrtemp,'Beam Coordinates Used'];
                case '10'
                    corstrtemp=[corstrtemp,'Instrument Coordinates Used'];
                case '01'
                    corstrtemp=[corstrtemp,'Ship Coordinates Used'];
                case '11'
                    corstrtemp=[corstrtemp,'Earth Coordinates Used'];
            end
            if corid(1)=='1' 
                corstrtemp=[corstrtemp,', Bin mapping used'];
            else
                corstrtemp=[corstrtemp,', Bin mapping not used'];                  
            end
            if corid(2)=='1'
                corstrtemp=[corstrtemp,', 3 beams solution allowed'];
            else
                corstrtemp=[corstrtemp,', 3 beams solution not allowed'];                  
            end
            if corid(3)=='1'
                corstrtemp=[corstrtemp,', Tilts used (only for ship or eart coordinates)'];
            else
                corstrtemp=[corstrtemp,', Tilts not used'];                  
            end
corstr=corstrtemp;

%%          FUNCTION TO CHECK DIFFERENCES IN FIXED LEADER
function issame=CheckFL(nValidFiles)
global dataout
issame=true;                                                               % Initialize result as true
if nValidFiles<2                                                           % If less than two files
    issame=true;                                                           % Return true
    return                                                                 % Exit function (Nothing to check!!)
end
for cntfiles=2:nValidFiles                                                 % Loop for all files
    issame=issame &&...                                                    % Check if any of the variables changes
        isequal(dataout.firmver(1),dataout.firmver(cntfiles))&&...
    isequal(dataout.firmrev(1),dataout.firmrev(cntfiles))&&...
    isequal(dataout.sysconf(1,:),dataout.sysconf(cntfiles,:))&&...
    isequal(dataout.SymData(1),dataout.SymData(cntfiles))&&...
    isequal(dataout.LagLength(1),dataout.LagLength(cntfiles))&&...
    isequal(dataout.usedbeams(1),dataout.usedbeams(cntfiles))&&...
    isequal(dataout.nbins(1),dataout.nbins(cntfiles))&&...
    isequal(dataout.pingperens(1),dataout.pingperens(cntfiles))&&...
    isequal(dataout.binsize(1),dataout.binsize(cntfiles))&&...
    isequal(dataout.blnk(1),dataout.blnk(cntfiles))&&...
    isequal(dataout.minthrsh(1),dataout.minthrsh(cntfiles))&&...
    isequal(dataout.ncodrep(1),dataout.ncodrep(cntfiles))&&...
    isequal(dataout.minpercgood(1),dataout.minpercgood(cntfiles))&&...
    isequal(dataout.maxerrvel(1),dataout.maxerrvel(cntfiles))&&...
    isequal(dataout.Tbetweenpng(1),dataout.Tbetweenpng(cntfiles))&&...
    isequal(dataout.corinfo(1,:),dataout.corinfo(cntfiles,:))&&...
    isequal(dataout.headalign(1),dataout.headalign(cntfiles))&&...
    isequal(dataout.headbias(1),dataout.headbias(cntfiles))&&...
    isequal(dataout.sensource(1,:),dataout.sensource(cntfiles,:))&&...
    isequal(dataout.senavail(1,:),dataout.senavail(cntfiles,:))&&...
    isequal(dataout.distmidbin1(1),dataout.distmidbin1(cntfiles))&&...
    isequal(dataout.lngthtranspulse(1),dataout.lngthtranspulse(cntfiles))&&...
    isequal(dataout.watrefbins(1,:),dataout.watrefbins(cntfiles,:))&&...
    isequal(dataout.mintarget(1),dataout.mintarget(cntfiles))&&...
    isequal(dataout.lowlattrig(1),dataout.lowlattrig(cntfiles))&&...
    isequal(dataout.distpulse(1),dataout.distpulse(cntfiles))&&...
    isequal(dataout.cpuserial(1,:),dataout.cpuserial(cntfiles,:))&&...
    isequal(dataout.bandwidth(1),dataout.bandwidth(cntfiles))&&...
    isequal(dataout.syspower(1),dataout.syspower(cntfiles))&&...
    isequal(dataout.basefreqid(1),dataout.basefreqid(cntfiles))&&...
    isequal(dataout.serial(1,:),dataout.serial(cntfiles,:))&&...
    isequal(dataout.HADCPbeamangle(1),dataout.HADCPbeamangle(cntfiles));
    if issame==false
        break
    end
end
