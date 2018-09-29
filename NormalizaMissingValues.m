function mDataN = NormalizaMissingValues(mData,media,desvia,R)

mDataN = mData;

for i = 1:size(mData,2)
    mDataN(R(:,i)~=1,i) = (mData(R(:,i)~=1,i) - media(i))/desvia(i);
end