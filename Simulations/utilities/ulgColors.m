classdef ulgColors
    %ULGCOLORS regroups all ulg colors as enum, the color can then be
    %expressed as a rgb array or a hex string
    properties
        R uint8
        G uint8
        B uint8
    end
    methods
        function obj = ulgColors(r, g, b)
            obj.R = r; obj.G = g; obj.B = b;
        end
    end
    enumeration
        % --- Primary colors
        BlueGreen (0,112,127)
        BlueGreenLight (95,164,176)
        BeigeLight (232,226,222)
        BeigePale (230,230,225)
        % --- Secondary colors (faculty)
        Beige (198,192,180)
        Yellow (255,208,0)
        YellowOrange (248,170,0)
        Orange (240,127,60)
        Red (230,45,49)
        GreenPale (185,205,118)
        GreenLight (125,185,40)
        Green (40,155,56)
        GreenDark (0,132,59)
        BlueLight (31,186,219)
        Blue (0,92,169)
        BlueViolet (91,87,162)
        Lavander (141,166,214)
        Purple (168,88,158)
        Violet (91,37,125)
        Grey (140,139,130)
        GreyLight (181,180,169)
        % --- White color (non official)
        White (255,255,255)
    end
    methods
        function rgb_array = rgb(obj)
            rgb_array = double([obj.R, obj.G, obj.B])/255;
        end
        function hex_string = hex(obj)
            hex_string = sprintf("#%02X%02X%02X",obj.R, obj.G, obj.B);
        end
    end
end

