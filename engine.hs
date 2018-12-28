module Engine where

{-# LANGUAGE BangPatterns #-}

import qualified Data.Map.Strict as Map
import Data.Map.Strict (Map)

import Data.List 

--(intercalate, sortOn, elemIndex)

import Data.Maybe (fromJust)

import System.Directory (getDirectoryContents)

import System.Random (randomRIO)

import System.IO

import Control.DeepSeq
import Control.Exception (evaluate)

-- ######## algebraic functions ########

type Position = (Float, Float, Float)

dot :: Position -> Position -> Float
dot (a1, a2, a3) (b1, b2, b3) = a1 * b1 + a2 * b2 + a3 * b3

reversedot :: Position -> Position -> Float -> Position
reversedot (a1, a2, a3) (b1, b2, b3) f = (a1 + b1 * cos f, a2 + b2 * cos f, a3 + b3 * cos f)

cross :: Position -> Position -> Position
cross (a1, a2, a3) (b1, b2, b3) = (a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1)

add :: Position -> Position -> Position
add (a1, a2, a3) (b1, b2, b3) = (a1 + b1, a2 + b2, a3 + b3)

sub :: Position -> Position -> Position
sub (a1, a2, a3) (b1, b2, b3) = (a1 - b1, a2 - b2, a3 - b3)

norm :: (Position, Position, Position) -> Position
norm (n, ca, c) = cross (sub ca n) (sub ca c)

dist :: Position -> Float
dist (a1, a2, a3) = sqrt $ a1^2 + a2^2 + a3^2

dihedral :: (Position, Position, Position) -> (Position, Position, Position) -> Float
dihedral a b = acos $ (dot v1 v2) / (dist v1 * dist v2)
    where v1 = norm a
          v2 = norm b 

unit :: Position -> Position
unit (a, b, c) = let s = a + b + c in (a / s, b / s, c / s)

-- ######## parsing ########

type Protein = [(String, (Float, Float))]

three2one :: String -> Char
three2one "ARG" = 'R'
three2one "HIS" = 'H'
three2one "LYS" = 'K'
three2one "ASP" = 'D'
three2one "GLU" = 'E'
three2one "SER" = 'S'
three2one "THR" = 'T'
three2one "ASN" = 'N'
three2one "GLN" = 'Q'
three2one "CYS" = 'C'
three2one "GLY" = 'G'
three2one "PRO" = 'P'
three2one "ALA" = 'A'
three2one "VAL" = 'V'
three2one "ILE" = 'I'
three2one "LEU" = 'L'
three2one "MET" = 'M'
three2one "PHE" = 'F'
three2one "TYR" = 'Y'
three2one "TRP" = 'W'
three2one _ = 'Z'

cif2prot :: String -> Protein
cif2prot s = toProt . collectAminos . filter (\x -> (not . null) x && head x == "ATOM") $ map words (lines s)
    where collectAminos (n:ca:c:xs) = 
            let 
                r = \(a, b, c) -> (read a :: Float, read b :: Float, read c :: Float)
                posN = r (n !! 10, n !! 11, n !! 12)
                posCA = r (ca !! 10, ca !! 11, ca !! 12)
                posC = r (c !! 10, c !! 11, c !! 12)
            in
                if three2one (n !! 5) == 'Z' then []
                else if n !! 3 == "N" && ca !! 3 == "CA" && c !! 3 == "C" then (n !! 5, (posN, posCA, posC)):collectAminos xs 
                else collectAminos (ca:c:xs)
          collectAminos _ = []
          toProt l | (not . null) l = ((fst . head) l, (0, 0)):toProtAux (tail l)
                   | otherwise = []
          toProtAux (a1:a2:a3:as) = (fst a2, (dihedral (snd a1) (snd a2), dihedral (snd a2) (snd a3))):toProtAux (a2:a3:as)
          toProtAux    (a2:a3:as) = [(fst a3, (0, 0))]
          toProtAux             _ = []

buildMap :: FilePath -> IO (Map String [((Float, Float), Int)])
buildMap fp =
    let
        merge l = Map.toList $ foldr (\x acc -> Map.insertWith (+) x 1 acc) Map.empty l
        scan (a1:a2:a3:as) = Map.insertWith (++) (map (three2one . fst) [a1, a2, a3]) [(snd a2)] $ scan (a2:a3:as)
        scan _ = Map.empty
        process [] = return Map.empty
        process (x:xs) =
            do
                p <- cif2prot <$> readFile x
                evaluate (force p)
                Map.unionWith (++) (scan p) <$> process xs
    in
        do
            dir <- filter (\x -> length x == 2) <$> getDirectoryContents fp
            l <- concat <$> mapM (\x -> map (\f -> (fp ++ x ++ "/" ++ f)) . filter (\f -> length f == 8) <$> getDirectoryContents (fp ++ x)) dir
            let l' = take 100000 l
            Map.map merge <$> process l'

avg' :: [((Float, Float), Int)] -> (Float, Float)
avg' l = (s1 / s3, s2 / s3)
    where s1 = sum $ map (\x -> fst (fst x) * fromIntegral (snd x)) l
          s2 = sum $ map (\x -> snd (fst x) * fromIntegral (snd x)) l
          s3 = fromIntegral . sum $ map snd l

var' :: [((Float, Float), Int)] -> Float
var' l = let a = avg' l in (foldr (\((a1, a2), _) acc -> acc + (a1 + a2 - fst a - snd a)^2) 0 l) / fromIntegral (length l)

closestN :: Int -> Int -> [((Float, Float), Int)] -> [((Float, Float), Int)]
closestN n i l = let e = fst (l !! i) in take n $ sortOn (\x -> sqrt $ (fst (fst x) - fst e)^2 + (snd (fst x) - snd e)^2) l

nclusters :: Int -> Map String [((Float, Float), Int)] -> Conditionals
nclusters n m = Map.map (map cast . repeat n . sort) m
    where sort = sortOn snd
          cast = \(a, b) -> ((round ((180 / pi) * fst a), round ((180 / pi) * snd a)), b)
          repeat n l | length l > n = let c = clump 10 l in repeat n $ (avg' c, sum (map snd c)):(l \\ c)
                     | otherwise = l
          clump n l = (\i -> closestN n i l) . (fst . head) . sortOn (snd) . snd $
                        foldr (\x (i, acc) -> (i+1, (i,x):acc)) (0, []) $ 
                            map (\x -> var' $ closestN n (fromJust (elemIndex x l)) l) l

-- model

toInt :: Protein -> [(String, (Int, Int))]
toInt l = map (\(s, (a, b)) -> (s, (round ((180 / pi) * a), round ((180 / pi) * b)))) l

toFloat :: [(String, (Int, Int))] -> Protein
toFloat l = map (\(s, (a, b)) -> (s, ((pi / 180) * fromIntegral a, (pi / 180) * fromIntegral b))) l

type Conditionals = Map String [((Int, Int), Int)]

loss :: Protein -> Protein -> Float
loss [] _ = 0
loss _ [] = 0
loss ((_, (a1, a2)):ps) ((_, (b1, b2)):qs) = loss ps qs + (a1 - b1)^2 + (a2 - b2)^2

avg :: [Float] -> Float
avg l = sum l / fromIntegral (length l)

var :: [Float] -> Float
var l = (foldr (\x acc -> acc + (x - a)^2) 0 l) / fromIntegral (length l)
    where a = avg l

sample :: [a] -> IO a
sample l = (\n -> (l !! n)) <$> randomRIO (0, length l - 1)

score :: Protein -> Float
score p = foldr (\x acc -> acc + foldr (\y acc' -> if dist (sub x y) < 5 then acc' + 1 else acc') 0 m) 0 m
    where m = Map.keys . (\(_, _, _, c) -> c) $ foldr (\x (pos, v, n, map) -> let (a, b, c) = update pos (snd x) v n in (a, b, c, Map.insert pos (fst x) map)) ((0, 0, 0), (1, 0, 0), (0, 0, 1), Map.empty) p
          update p (a, b) v n = (add p v, unit (cross n v), reversedot n v a)

--score p = foldr (\x (pos, map) -> Map.insert (fst x) () map) ((0, 0, 0), Map.empty) p

generateMLE :: Int -> Protein -> Conditionals -> IO [(String, (Int, Int))]
generateMLE 5 p m | length p > 2 = (\x -> ((fst . head) p, (0, 0)):((fst . head . tail) p, (0, 0)):x) <$> generateMLE5 p m
generateMLE 7 p m | length p > 3 = (\x -> ((fst . head) p, (0, 0)):((fst . head . tail) p, (0, 0)):((fst . head . tail . tail) p, (0, 0)):x) <$> generateMLE7 p m
generateMLE _ p m | length p > 0 = (((fst . head) p, (0, 0)):) <$> generateMLE3 p m
generateMLE _ _ _ = return []

generateMLE3 :: Protein -> Conditionals -> IO [(String, (Int, Int))]
generateMLE3 (a1:a2:a3:as) m = 
    let
        str = map (three2one . fst) [a1, a2, a3]
        ang = fst . foldr (\x acc -> if snd x > snd acc then x else acc) ((0, 0), 0) $ fromJust (Map.lookup str m)
    in
        do
            ((fst a2, ang):) <$> generateMLE3 (a2:a3:as) m 
generateMLE3 (a2:a3:[]) m = return [(fst a3, (0, 0))]
generateMLE3 _ _ = return []

generateMLE5 :: Protein -> Conditionals -> IO [(String, (Int, Int))]
generateMLE5 (a1:a2:a3:a4:a5:as) m = 
    let
        str3 = map (three2one . fst) [a2, a3, a4]
        str5 = map (three2one . fst) [a1, a2, a3, a4, a5]
        ang3 = fst . foldr (\x acc -> if snd x > snd acc then x else acc) ((0, 0), 0) <$> Map.lookup str3 m
        ang5 = fst . foldr (\x acc -> if snd x > snd acc then x else acc) ((0, 0), 0) <$> Map.lookup str5 m
    in
        do
            ang <- case ang5 of
                       Just a5 -> return a5
                       Nothing -> case ang3 of
                                      Just a3 -> return a3
                                      Nothing -> return (0, 0)
            ((fst a3, ang):) <$> generateMLE5 (a2:a3:a4:a5:as) m 
generateMLE5 (a2:a3:a4:a5:[]) m = generateMLE3 [a3, a4, a5] m
generateMLE5 _ _ = return []

generateMLE7 :: Protein -> Conditionals -> IO [(String, (Int, Int))]
generateMLE7 (a1:a2:a3:a4:a5:a6:a7:as) m = 
    let
        str3 = map (three2one . fst) [a3, a4, a5]
        str5 = map (three2one . fst) [a2, a3, a4, a5, a6]
        str7 = map (three2one . fst) [a1, a2, a3, a4, a5, a6, a7]
        ang3 = fst . foldr (\x acc -> if snd x > snd acc then x else acc) ((0, 0), 0) <$> Map.lookup str3 m
        ang5 = fst . foldr (\x acc -> if snd x > snd acc then x else acc) ((0, 0), 0) <$> Map.lookup str5 m
        ang7 = fst . foldr (\x acc -> if snd x > snd acc then x else acc) ((0, 0), 0) <$> Map.lookup str7 m
    in
        do
            ang <- case ang7 of
                       Just a7 -> return a7
                       Nothing -> case ang5 of
                                      Just a5 -> return a5
                                      Nothing -> case ang3 of
                                                     Just a3 -> return a3
                                                     Nothing -> return (0, 0)
            ((fst a4, ang):) <$> generateMLE7 (a2:a3:a4:a5:a6:a7:as) m 
generateMLE7 (a2:a3:a4:a5:a6:a7:[]) m = generateMLE5 [a3, a4, a5, a6, a7] m
generateMLE7 _ _ = return []

generate :: Int -> Protein -> Conditionals -> IO [(String, (Int, Int))]
generate 5 p m | length p > 2 = (\x -> ((fst . head) p, (0, 0)):((fst . head . tail) p, (0, 0)):x) <$> generate5 p m
generate 7 p m | length p > 3 = (\x -> ((fst . head) p, (0, 0)):((fst . head . tail) p, (0, 0)):((fst . head . tail . tail) p, (0, 0)):x) <$> generate7 p m
generate _ p m | length p > 0 = (((fst . head) p, (0, 0)):) <$> generate3 p m
generate _ _ _ = return []


generate3 :: Protein -> Conditionals -> IO [(String, (Int, Int))]
generate3 (a1:a2:a3:as) m = 
    let
        str = map (three2one . fst) [a1, a2, a3]
    in
        do
            ang <- case Map.lookup str m of
                       Just l -> sample $ concatMap (\x -> replicate (snd x) (fst x)) l
                       Nothing -> return (0, 0)
            ((fst a2, ang):) <$> generate3 (a2:a3:as) m 
generate3 (a2:a3:[]) m = return [(fst a3, (0, 0))]
generate3 _ _ = return []

generate5 :: Protein -> Conditionals -> IO [(String, (Int, Int))]
generate5 (a1:a2:a3:a4:a5:as) m = 
    let
        str3 = map (three2one . fst) [a2, a3, a4]
        str5 = map (three2one . fst) [a1, a2, a3, a4, a5]
        l3 = Map.lookup str3 m
        l5 = Map.lookup str5 m
        --ang3 = sample $ fromJust . concatMap (\x -> replicate (snd x) (fst x)) <$> l3
        --ang5 = sample $ fromJust . concatMap (\x -> replicate (snd x) (fst x)) <$> l5
        --score3 = (\l -> (sum . (map snd)) l) . fromJust <$> l3
        --score5 = sum . (map snd) <$> Map.lookup str5 m
    in
        do
            ang <- case Map.lookup str5 m of
                       Just l -> sample $ concatMap (\x -> replicate (snd x) (fst x)) l
                       Nothing -> case Map.lookup str3 m of
                                      Just l' -> sample $ concatMap (\x -> replicate (snd x) (fst x)) l'
                                      Nothing -> return (0, 0)
            ((fst a3, ang):) <$> generate5 (a2:a3:a4:a5:as) m 
generate5 (a2:a3:a4:a5:[]) m = generate3 [a3, a4, a5] m
generate5 _ _ = return []

generate7 :: Protein -> Conditionals -> IO [(String, (Int, Int))]
generate7 (a1:a2:a3:a4:a5:a6:a7:as) m = 
    let
        str3 = map (three2one . fst) [a3, a4, a5]
        str5 = map (three2one . fst) [a2, a3, a4, a5, a6]
        str7 = map (three2one . fst) [a1, a2, a3, a4, a5, a6, a7]
    in
        do
            ang <- case Map.lookup str7 m of
                       Just l -> sample $ concatMap (\x -> replicate (snd x) (fst x)) l
                       Nothing -> case Map.lookup str5 m of
                                      Just l' -> sample $ concatMap (\x -> replicate (snd x) (fst x)) l'
                                      Nothing -> case Map.lookup str3 m of
                                                 Just l'' -> sample $ concatMap (\x -> replicate (snd x) (fst x)) l''
                                                 Nothing -> return (0, 0)
            ((fst a4, ang):) <$> generate7 (a2:a3:a4:a5:a6:a7:as) m 
generate7 (a2:a3:a4:a5:a6:a7:[]) m = generate5 [a3, a4, a5, a6, a7] m
generate7 _ _ = return []

-- testing

prior :: IO (Map String [((Int, Int), Int)])
prior = read <$> readFile "comp_dist_big.txt"


{-markovchain :: Protein -> Conditionals -> IO [(String, (Int, Int))]
markovchain p m = 
    let
        replace i x l = take i l ++ x:drop (i+1) l
        sampleProt s = 
            do
                i <- randomRIO (0, length s - 1)
                case Map.lookup (fst (s !! i)) m of
                    Nothing -> return s
                    Just l -> (\x -> replace (elemIndex x s) (fst (s !! i), x) s) <$> sample l
        repeat 0 s = return []
        repeat n s = 
            do
                s' <- toFloat <$> sampleProt s
                (\x -> (s, score (toFloat s)):x) <$> repeat (n-1) s'
    in
        do 
            state <- generateMLE7 p m 
            l <- repeat 100 state
            return fst $ foldr (\x acc -> if snd x > snd acc then x else acc) (head l) l 
-}

main :: IO ()
main = 
    let
        repeat _ 0 _ _ = return []
        repeat f n p m =
            do
                p' <- toFloat <$> f p m
                evaluate (force p')
                (((loss p p') / fromIntegral (length p)):) <$> repeat f (n-1) p m
        process'' _ _ [] _ = return []
        process'' f n (x:xs) m = 
            do 
                p <- cif2prot <$> readFile x
                evaluate (force p)
                p' <- repeat f n p m 
                if length p < 3 then process'' f n xs m
                else ((minimum p'):) <$> process'' f n xs m
        process' _ _ [] _ = return []
        process' f n (x:xs) m = 
            do 
                p <- cif2prot <$> readFile x
                evaluate (force p)
                p' <- repeat f n p m 
                if length p < 3 then process' f n xs m
                else ((avg p'):) <$> process' f n xs m
        process _ _ [] _ = return []
        process f n (x:xs) m = 
            do 
                p <- cif2prot <$> readFile x
                evaluate (force p)
                p' <- toFloat <$> f p m 
                evaluate (force p')
                if length p < n then process f n xs m
                else (((loss p p') / fromIntegral (length p)):) <$> process f n xs m
    in 
        do
            m <- prior
            evaluate (force m)
            dir <- filter (\x -> length x == 2) <$> getDirectoryContents "mmCIF/"
            l <- concat <$> mapM (\x -> map (\f -> ("mmCIF/" ++ x ++ "/" ++ f)) . filter (\f -> length f == 8) <$> getDirectoryContents ("mmCIF/" ++ x)) dir
            let l' = take 500 l
            
            error3 <- process (generate3) 3 l' m
            error5 <- process (generate5) 5 l' m
            error7 <- process (generate7) 7 l' m
            print $ "Error 3: Avg. " ++ show (avg error3) ++ " Var. " ++ show (var error3)
            print $ "Error 5: Avg. " ++ show (avg error5) ++ " Var. " ++ show (var error5)
            print $ "Error 7: Avg. " ++ show (avg error7) ++ " Var. " ++ show (var error7)

            error30 <- process' (generate3) 10 l' m
            error50 <- process' (generate5) 10 l' m
            error70 <- process' (generate7) 10 l' m
            error30' <- process'' (generate3) 10 l' m
            error50' <- process'' (generate5) 10 l' m
            error70' <- process'' (generate7) 10 l' m
            print $ "Error 30: Avg. " ++ show (avg error30) ++ " Var. " ++ show (var error30)
            print $ "Error 50: Avg. " ++ show (avg error50) ++ " Var. " ++ show (var error50)
            print $ "Error 70: Avg. " ++ show (avg error70) ++ " Var. " ++ show (var error70)
            print $ "Error best of 30: Avg. " ++ show (avg error30') ++ " Var. " ++ show (var error30')
            print $ "Error best of 50: Avg. " ++ show (avg error50') ++ " Var. " ++ show (var error50')
            print $ "Error best of 70: Avg. " ++ show (avg error70') ++ " Var. " ++ show (var error70')

            error300 <- process' (generate3) 100 l' m
            error500 <- process' (generate5) 100 l' m
            error700 <- process' (generate7) 100 l' m
            error300' <- process'' (generate3) 100 l' m
            error500' <- process'' (generate5) 100 l' m
            error700' <- process'' (generate7) 100 l' m
            print $ "Error 300: Avg. " ++ show (avg error300) ++ " Var. " ++ show (var error300)
            print $ "Error 500: Avg. " ++ show (avg error500) ++ " Var. " ++ show (var error500)
            print $ "Error 700: Avg. " ++ show (avg error700) ++ " Var. " ++ show (var error700)
            print $ "Error best of 300: Avg. " ++ show (avg error300') ++ " Var. " ++ show (var error300')
            print $ "Error best of 500: Avg. " ++ show (avg error500') ++ " Var. " ++ show (var error500')
            print $ "Error best of 700: Avg. " ++ show (avg error700') ++ " Var. " ++ show (var error700')
            
            mle3 <- process (generateMLE3) 3 l' m
            mle5 <- process (generateMLE5) 5 l' m
            mle7 <- process (generateMLE7) 7 l' m
            print $ "MLE Error 3: Avg. " ++ show (avg mle3) ++ " Var. " ++ show (var mle3)
            print $ "MLE Error 5: Avg. " ++ show (avg mle5) ++ " Var. " ++ show (var mle5)
            print $ "MLe Error 7: Avg. " ++ show (avg mle7) ++ " Var. " ++ show (var mle7)
            





