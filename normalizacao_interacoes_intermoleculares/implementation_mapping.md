# Mapeamento implementado no DockLens

## Perfis científicos

| Chave interna | Interface | Política |
|---|---|---|
| `plip` | PLIP (legacy default) | Preserva o comportamento histórico. |
| `luna` | LUNA 0.14 defaults | Usa os limites nativos registrados em `luna_defaults.tsv` e exige hidrogênio explícito para doadores. |
| `dsv` | DSV-like (Discovery Studio 2024 defaults) | Usa os critérios retornados por `Mdm::NonbondMonitor` e percepção química conservadora. |
| `luna_dsv` | LUNA × DSV conservative | Une famílias representáveis e aplica o critério compartilhado mais restritivo. |

O nome persistido `hbond_preset` foi mantido para compatibilidade com projetos e manifests anteriores, mas agora identifica o perfil científico completo.

## Regra do perfil cruzado

- distância máxima compartilhada: menor valor;
- ângulo mínimo compartilhado: maior valor;
- ângulo máximo compartilhado: menor valor;
- geometrias adicionais de uma das fontes permanecem obrigatórias quando aplicáveis;
- halogênio e calcogênio precisam satisfazer simultaneamente distância absoluta e fração VDW quando os dois modelos fornecem essas grandezas;
- pares carga oposta dentro do limite de ponte salina são classificados somente como `saltbridge`; `attractive_charge` cobre apenas a faixa adicional, evitando duplicação;
- contatos atômicos genéricos `proximal`, `vdw` e `vdw_clash` não são habilitados por padrão;
- famílias com múltiplas combinações atômicas são reduzidas ao contato representativo mais próximo por par semântico de resíduos/grupos.

## Famílias adicionais representadas

- `chalcogen`: modelo conservador R–Y···A–N para Y = S, Se ou Te;
- `attractive_charge`: atração eletrostática fora do núcleo já classificado como ponte salina;
- `charge_repulsion`: contato entre centros de carga de mesmo sinal.

## Autoridade de cores

`docklens.interaction_core.INTERACTION_COLORS` é a autoridade canônica. Tabelas, gráficos, planilhas e o plugin PyMOL consomem a mesma ordem, rótulos e cores. Como há mais famílias que cores qualitativas Okabe–Ito, o plugin preserva a cor semântica e usa padrões pontilhados para as famílias estendidas.

## Limites da equivalência

O perfil LUNA porta os parâmetros públicos para a percepção química sem dependências do DockLens; ele não incorpora Open Babel nem promete saída idêntica ao executável LUNA. O perfil DSV-like reproduz os critérios expostos pela API instalada, não algoritmos proprietários não publicados.
