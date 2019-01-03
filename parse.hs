module HPModel where

import Control.Applicative

import qualified Data.Map.Strict as Map
import Data.Map.Strict (Map)

import qualified Data.Set as Set
import Data.Set (Set)

import qualified Data.Char as Char
import Data.Char (Char)

import qualified System.Random as Random
import System.Random (Random)

-- data

data Atom = Bridge Int
          | AlphaHelix String
          | BetaStrand String 
          | BetaSheet String 
          | Chain Int
          | Glob String

type BondGraph = Map Int [Int]

instance Show (Atom) where
    show (Bridge n) = "Bridge " ++ show n
    show (Chain n) = "Chain " ++ show n
    show (AlphaHelix s) = "AlphaHelix " ++ show (length s)
    show (BetaSheet s) = "BetaSheet " ++ show (length s)
    show (Glob s) = "Glob " ++ show (length s) ++ " " ++ show (length (filter polar s))

newtype Parser a = Parser {parse :: String -> Maybe (a, String)}

instance Functor Parser where
    fmap f p = Parser $ \s -> (\(a, c) -> (f a, c)) <$> parse p s

instance Applicative Parser where
    pure a = Parser $ \s -> Just (a, s)
    f <*> a = Parser $ \s -> 
        case parse f s of
            Just (g, s') -> parse (fmap g a) s'
            Nothing -> Nothing

instance Alternative Parser where
    empty = Parser $ \_ -> Nothing
    l <|> r = Parser $ \s -> parse l s <|> parse r s 

-- helper functions

amino :: Char -> Bool
amino 'r' = True
amino 'h' = True
amino 'k' = True
amino 'd' = True
amino 'e' = True
amino 's' = True
amino 't' = True
amino 'n' = True
amino 'q' = True
amino 'c' = True
amino 'g' = True
amino 'p' = True
amino 'a' = True
amino 'v' = True
amino 'i' = True
amino 'l' = True
amino 'm' = True
amino 'f' = True
amino 'y' = True
amino 'w' = True
amino _   = False

charge :: Char -> Int
charge 'r' = 1
charge 'h' = 1
charge 'k' = 1
charge 'd' = -1
charge 'e' = -1
charge _   = 0

polar :: Char -> Bool
polar 's' = True
polar 't' = True
polar 'n' = True
polar 'q' = True
polar c   = charge c /= 0

count :: (a -> Bool) -> [a] -> Int
count f (x:xs) | f x = 1 + count f xs
count _ _ = 0

ensure :: (a -> Bool) -> Parser a -> Parser a
ensure p parser = Parser $ \s ->
    case parse parser s of
        Nothing -> Nothing
        Just (a,s') -> if p a then Just (a,s') else Nothing

-- parsers

bridge :: Parser Int
bridge = Parser $ \s -> 
    case count (=='c') s of 
        0 -> Nothing
        n -> Just (n, drop n s)

chain :: Parser Int
chain = Parser $ \s -> 
    case count (\c -> c =='p' || c == 'g') s of 
        0 -> Nothing
        n -> Just (n, drop n s)

alphahelix :: Parser String
alphahelix = Parser $ \s -> let n = count polar s
                            in if n > 5 then Just (take n s, drop n s) else Nothing 

betaturn :: Parser ()
betaturn = Parser f 
    where f ('p':'p':xs) = Just ((), xs)
          f ('g':'g':xs) = Just ((), xs)
          f ('p':'g':xs) = Just ((), xs)
          f ('g':'p':xs) = Just ((), xs)
          f _ = Nothing

betastrand :: Parser Int
betastrand = let p = Parser f where f (x:y:xs) | polar x /= polar y = Just (1, y:xs)
                                    f _ = Nothing  
             in Parser $ \s ->
                case parse (some p) s of 
                    Nothing -> Nothing
                    Just (n, s) -> if (sum n) < 3 then Nothing else Just (sum n, s)
    
betasheet :: Parser String
betasheet = Parser $ \s -> 
    case parse (some p) s of
        Nothing -> Nothing
        Just (l, s') -> if f l (head l) && sum l > 9 then Just (take (sum l) s, s') else Nothing
    where p = betastrand <* betaturn 
          f lst n = foldr (\x acc -> acc && x == n) True lst

glob :: Parser String
glob = Parser $ \s -> 
    let 
        collect _ [] = []
        collect p (x:xs) = case parse p (x:xs) of
                                Nothing -> x:collect atom xs
                                Just _ -> [] 
        l = collect atom s
    in
        if l /= [] then Just (l, drop (length l) s)  else Nothing

atom :: Parser Atom
atom = Bridge <$> bridge
   <|> AlphaHelix <$> alphahelix
   <|> BetaSheet <$> betasheet
   <|> Chain <$> chain

term :: Parser Atom
term = Glob <$> glob 
   <|> atom

backbone :: Parser [Atom]
backbone = some term

protein :: Parser BondGraph
protein = undefined

type Partition = (Set Int, Set Int)

{--collect :: Int -> Backbone -> Partition
collect l = 
    let 
        f i lst = foldl ((j, h, p) x -> if polar x then (j + 1, (h, Set.insert j p)) else (j + 1, (Set.insert i h, p))) (i, (Set.empty, Set.empty)) lst
        help n [] = Set.empty
        help n (x:xs) = case x of
                            Glob g -> Set.union (f n g) (help (n + length g + 1) xs)
                            Bridge m -> help (n + m + 1) xs
                            Chain m -> help (n + m + 1) xs
                            AlphaHelix s -> help (n + length s + 1) xs
                            BetaSheet s -> help (n + length s + 1) xs 
    in
        help 0 l

collect :: Int -> Backbone -> Set Int
collect n [] = Set.empty
collect n (x:xs) = case x of
                        Glob g -> Set.union (Set.fromList [n..(n + length g)]) (collect (n + length g + 1) xs)
                        Bridge m -> collect (n + m) xs
                        Chain m -> collect (n + m) xs
                        AlphaHelix s -> collect (n + length s) xs
                        BetaSheet s -> collect (n + length s) xs --}

-- HP Model

data HP = H | P | G | C deriving (Eq, Show)

type LookupTable = Map Int HP

type ContactMap = Map Int [Int]

y :: String -> Int
y [] = 0
y ('y':xs) = 1 + y xs
y (_:xs) = y xs

f :: String -> Int
f [] = 0
f ('f':xs) = 1 + f xs
f (_:xs) = f xs

big :: String -> Int
big [] = 0
big ('w':xs) = 1 + big xs
big ('y':xs) = 1 + big xs
big ('f':xs) = 1 + big xs
big (_:xs) = big xs

len :: Char -> Int
len c | (not . amino ) c = 0
len 'w' = 2
len 'y' = 2
len 'f' = 2
len 'm' = 2
len 'l' = 1
len 'i' = 1
len 'v' = 1
len 'a' = 1
len 'k' = 2
len 'h' = 2
len 'r' = 2

lexer :: String -> [HP] 
lexer [] = []
lexer ('g':xs) = G:lexer xs
lexer ('c':xs) = C:lexer xs
lexer (x:xs) | (not . amino) x = lexer xs
             | polar x   = P:lexer xs
             | otherwise = H:lexer xs

reduce :: [HP] -> [(HP, Int)]
reduce [] = []
reduce s = if h /= 0 then (H, h):reduce (drop h s) 
      else if p /= 0 then (P, p):reduce (drop p s) 
      else if g /= 0 then (G, g):reduce (drop g s)
      else (C, c):reduce(drop c s)
    where h = count (==H) s
          p = count (==P) s
          g = count (==G) s
          c = count (==C) s

globule :: [(HP, Int)] -> [Int]
globule l = foldr (\x acc -> snd x:acc) [] $ filter (\(hp, n) -> hp == H) l

helix :: [(HP, Int)] -> [Int]
helix [] = []
helix ((P, n):xs) | n > 2 = n:helix xs
helix (_:xs) = helix xs


h :: Int -> [(HP, Int)] -> [Int]
h _ [] = []
h n ((H, m):xs) = [n..(n + m - 1)] ++ (h (n + m) xs)
h n ((_, m):xs) = h (n + m) xs


data AminoAcid = 
    Arginine 
  | Histidine 
  | Lysine 
  | Aspartate 
  | Glutamate 
  | Serine 
  | Threonine 
  | Asparagine 
  | Glutamine 
  | Cysteine 
  | Glycine 
  | Proline
  | Alanine 
  | Valine
  | Isoleucine 
  | Leucine 
  | Methionine 
  | Phenylalanine 
  | Tyrosine 
  | Tryptophan

-- interpreter


globulin2u :: String
globulin2u =  "MKLLLLLLCLGLTLVCGHAEEASSTRGNLDVAKLNGDWFSIVVASNKREKBTTTTGGGGEEEEEEEESSGGGIEENGSMRVFMQHIDVLENSLGFKFRIKENGECRELYLVAYKTPEDGEYFTSTTTTEEEEEEEETTEEEEEEEEESSSEEEEEEEEEESSTTEEVEYDGGNTFTILKTDYDRYVMFHLINFKNGETFQLMVLYGRTKDLSSDIKEESSSEEEEEEEEESSSEEEEEEEEEETTEEEEEEEEEESSSSHHHH151EKFAKLCEAHGITRDNIIDLTKTDRCLQARGHHHHHHHHTTTGGGEEEGGGS"

p1nf3 :: String
p1nf3 = "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHHHHHHHHHHHGGGTGGGTTHHHHHHHHHHHHHSSSSSSSSSSTGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLTTTHHHHHHHHTHHHHHHHHGGGGGGGTSHHHHHHTTTTHHHHHTLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYRTTTTTTHHHHTTSTTTHHHHHHHHHHHHHHHHHSGGG"

s :: String
s = foldr (\x acc -> Char.toLower x:acc) [] globulin2u

s' :: String
s' = foldr (\x acc -> Char.toLower x:acc) [] p1nf3

tpi :: String
tpi = "MAEDGEEAEFHFAALYISGQWPRLRADTDLQRLGSSAMAPSRKFFVGGNWKMNGRKQSLGELIGTLNAAKVPADTEVVCAPPTAYIDFARQKLDPKIAVAAQNCYKVTNGAFTGEISPGMIKDCGATWVVLGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFEQTKVIADNVKDWSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLRGWLKSNVSDAVAQSTRIIYGGSVTGATCKELASQPDVDGFLVGGASLKPEFVDIINAKQ"

t :: String
t = map Char.toLower tpi






