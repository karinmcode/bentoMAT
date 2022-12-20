function fn0=my_find_wildcard_filename(url_with_wildcard)
%% fn0=my_find_wildcard_filename(url_with_wildcard)

if contains(url_with_wildcard,'*')
    files = dir(string(url_with_wildcard));%dir('Y:\Users\Karin\data\processed\aligned2vid\m943\220413\003\CAM1_mm943_220413_003*.mp4')
    % dir('Y:\Users\Karin\data\processed\aligned2vid\m532\211116\002\m532_211116_002*.mp4')
    fn0 = myfilename(files.name);
else
    fn0 = myfilename(url_with_wildcard);
end
