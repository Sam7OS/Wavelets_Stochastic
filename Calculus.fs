namespace Stocalc

module Calculus =
    open System
    
    let rnd = Random()

(*############################################################ - Types - ##################################################################*)

    type ParameterType = Constant of float | Deterministic of (float -> float) | Stochastic of (float*float) list
    type ProcessParameter = Drift of ParameterType | Diffusion of ParameterType 

    type Std_BM_params = { Time_Start:float ; Time_End:float ; Time_Scale:int }
    type Gen_BM_params = { Time_Start:float ; Time_End:float ; Time_Scale:int ; Start_Point:float ; Mu:ProcessParameter ; Sigma:ProcessParameter}
    type Geo_BM_params = { Time_Start:float ; Time_End:float ; Time_Scale:int ; Start_Point:float ; Mu:ProcessParameter ; Sigma:ProcessParameter}
    
    type Brownian = Std_Brownian of Std_BM_params | Gen_Brownian of Gen_BM_params | Geo_Brownian of Geo_BM_params

(*############################################################ - Probability Distributions - ##################################################################*)

    // Générateur de nombres aléatoires 
    let genUnifRandomNumbers length =
        List.init length (fun _ -> rnd.NextDouble())

    // Echantillon de Bernoulli
    let bernoulli_sample length p=
        genUnifRandomNumbers length |> List.map (fun v -> if v < p then 1 else 0)

    // Echantillon de loi binomiale
    let binomial_sample length n p=
        [1..n] |> List.map (fun _ -> (bernoulli_sample length p |> List.sum)) 

    // Echantillon de loi exponentielle
    let exponential_sample length intensity = 
        genUnifRandomNumbers length |> List.map (fun v -> -Math.Log(v)/intensity)

    // Echantillon de loi de Poisson
    let poisson_sample n intensity =
        let expo_sample = exponential_sample n intensity
        (expo_sample.Head,expo_sample)
        ||> List.scan (fun s x -> s + x)
        |> List.pairwise 

    // Distribution de Poisson
    let poisson_distribution intensity times =
        let poi_sample = poisson_sample (List.length times) intensity
        times 
        |> List.map (fun t -> poi_sample |> List.mapi (fun i z -> if (t < snd z && t > fst z) then i+1 else 0) |> List.sum)
        |> List.zip times

    // Echantillon de loi normale centrée réduite
    let gauss_centered_reduced length =
        let R = genUnifRandomNumbers length |> List.map (fun x -> Math.Sqrt(-2.0*Math.Log(x)))
        let Theta = genUnifRandomNumbers length |> List.map (fun x -> Math.Cos(2.0*Math.PI*x))
        List.map2 (fun x y -> x*y) R Theta

    // Echantillon de loi normale
    let gauss_list length mean std =
        gauss_centered_reduced length |> List.map (fun x -> mean + std*x )

    // Intégrale de la fonction anyFunc entre x0 et x1
    let Integral x0 x1 scale anyFunc = 
        (x1-x0) * (genUnifRandomNumbers scale |> List.map (fun v -> anyFunc (x0 + v * (x1-x0))) |> List.average)

    // Densité de loi normale
    let gauss_density mean std x =
        Math.Exp(-0.5*Math.Pow((x-mean)/std,2.0))/(std*Math.Sqrt(2.0*Math.PI))

    // Densité de loi normale centrée réduite
    let gauss_01_density x =
        gauss_density 0.0 1.0 x

    // Distribution de Poisson composé
    let compound_poisson_distribution intensity times =
        let poi = poisson_distribution intensity times
        let X,Y = poi |> List.unzip
        let last = Y |> List.last
        let states = List.init last (fun v -> v)
        let binomial = binomial_sample 5 last 0.4
        let acc_bino = (0,binomial) ||> List.scan (fun acc v -> acc + v)
        let res = Y |> List.map (fun v -> ( acc_bino |> List.mapi (fun i z -> if v = i then z else 0) ) |> List.sum )
        (X,res) ||> List.zip 

(*############################################################ - Stochastic Process - ##################################################################*)

    // Marche aléatoire
    let BrownianM_RdWalk time_start time_end time_scale start_point =
        let coin_res = start_point :: List.init time_scale (fun _ -> if rnd.Next(0,2) = 0 then -1.0 else 1.0) 
        let times = [0.0..float(time_scale)] |> List.map (fun v -> time_start + v*(time_end-time_start)/float(time_scale))
        List.zip times ((coin_res.Head,coin_res.Tail) ||> List.scan (fun acc x -> acc + x))

    // Mouvement brownien standard
    let StandardBrownianMotion time_start time_end time_scale = 
        let times = [0.0..float(time_scale)] |> List.map (fun v -> time_start + v*(time_end-time_start)/float(time_scale))
        (0.0 , gauss_list time_scale 0.0 1.0)
        ||> List.scan (fun acc v -> acc + Math.Sqrt((time_end-time_start)/float(time_scale))*v) 
        |> List.zip times

    // Mouvement Brownien
    let BrownianMotion time_start time_end time_scale start_point drift diffusion = 
        match drift , diffusion with
        // Constant Drift
        | Drift (Constant const_drift) , Diffusion (Constant const_diff) -> 
            Some(StandardBrownianMotion time_start time_end time_scale 
            |> List.map (fun v -> ( (fst v) , start_point + const_drift * (fst v) + const_diff * (snd v) )))
        | Drift (Constant const_drift) , Diffusion (Deterministic det_diff) ->
            Some(StandardBrownianMotion time_start time_end time_scale 
            |> List.map (fun v -> ( (fst v) , start_point + const_drift * (fst v) + det_diff(fst v) * (snd v) )))
        | Drift (Constant const_drift) , Diffusion (Stochastic stoch_diff) ->
            Some((stoch_diff , StandardBrownianMotion time_start time_end time_scale) 
            ||> List.map2 (fun z v -> ( (fst v) , start_point + const_drift * (fst v) + (snd z) * (snd v) )))
        // Deterministic Drift
        | Drift (Deterministic det_drift) , Diffusion (Constant const_diff) -> 
            Some(StandardBrownianMotion time_start time_end time_scale 
            |> List.map (fun v -> ( (fst v) , start_point + det_drift(fst v) * (fst v) + const_diff * (snd v) )))
        | Drift (Deterministic det_drift) , Diffusion (Deterministic det_diff) ->
            Some(StandardBrownianMotion time_start time_end time_scale 
            |> List.map (fun v -> ( (fst v) , start_point + det_drift(fst v) * (fst v) + det_diff(fst v) * (snd v) )))
        | Drift (Deterministic det_drift) , Diffusion (Stochastic stoch_diff) ->
            Some((stoch_diff , StandardBrownianMotion time_start time_end time_scale) 
            ||> List.map2 (fun z v -> ( (fst v) , start_point + det_drift(fst v) * (fst v) + (snd z) * (snd v) )))
        // Stochastic Drift
        | Drift (Stochastic stoch_drift) , Diffusion (Constant const_diff) -> 
            Some((stoch_drift , StandardBrownianMotion time_start time_end time_scale) 
            ||> List.map2 (fun z v -> ( (fst v) , start_point + (snd z) * (fst v) + const_diff * (snd v) )))
        | Drift (Stochastic stoch_drift) , Diffusion (Deterministic det_diff) ->
            Some((stoch_drift , StandardBrownianMotion time_start time_end time_scale) 
            ||> List.map2 (fun z v -> ( (fst v) , start_point + (snd z) * (fst v) + det_diff(fst v) * (snd v) )))
        | Drift (Stochastic stoch_drift) , Diffusion (Stochastic stoch_diff) ->
            Some((stoch_drift , stoch_diff , StandardBrownianMotion time_start time_end time_scale) 
            |||> List.map3 (fun z q v -> ( (fst v) , start_point + (snd z) * (fst v) + (snd q) * (snd v) )))    
        // Error
        | _ -> None

    // Mouvement Brownien géométrique
    let GeometricBrownianMotion time_start time_end time_scale start_point drift diffusion = 
        match drift , diffusion with
        // Constant Drift , Constant Diffusion
        | Drift (Constant const_drift) , Diffusion (Constant const_diff) ->
            Some(StandardBrownianMotion time_start time_end time_scale
            |> List.map (fun v -> ( (fst v) , start_point * Math.Exp((const_drift - 0.5 * Math.Pow(const_diff,2.0)) * (fst v) + const_diff * (snd v)) )))
        // Error
        | _ -> None

    // Permet la simulation d'un Mouvement Brownien
    let Create_BM brownian =
        match brownian with
        | Std_Brownian { Time_Start=time_start ; Time_End=time_end ; Time_Scale=time_scale } -> 
            if (time_start < time_end && time_start >= 0.0 && time_scale >= 1) then Some(StandardBrownianMotion time_start time_end time_scale) else None
        | Gen_Brownian { Time_Start=time_start ; Time_End=time_end ; Time_Scale=time_scale ; Start_Point=start_point ; Mu=drift ; Sigma=diffusion} -> 
            if (time_start < time_end && time_start >= 0.0 && time_scale >= 1) then BrownianMotion time_start time_end time_scale start_point drift diffusion else None
        | Geo_Brownian { Time_Start=time_start ; Time_End=time_end ; Time_Scale=time_scale ; Start_Point=start_point ; Mu=drift ; Sigma=diffusion} -> 
            if (time_start < time_end && time_start >= 0.0 && time_scale >= 1) then GeometricBrownianMotion time_start time_end time_scale start_point drift diffusion else None
        
       
(*############################################################ - CSV Creation - ##################################################################*)

    // Pour écrire la simulation dans un CSV avec la variable temps
    let write_csv_2col file_path list =
        match list with
        | Some list ->
            let csv = new System.IO.StreamWriter(path=file_path)
            list |> List.map (fun v -> csv.WriteLine(fst(v).ToString() + ";" + snd(v).ToString())) |> ignore
            csv.Close()
        | None -> None |> ignore

    // Pour écrire la simulation dans un CSV
    let write_csv_1col file_path list =
        match list with
        | Some list ->
            let csv = new System.IO.StreamWriter(path=file_path)
            list |> List.map (fun v -> csv.WriteLine(v.ToString()+";")) |> ignore
            csv.Close()
        | None -> None |> ignore