#!/usr/bin/env python3
"""Tool groups for the retained interactive research runtime.

Groups contain only the built-in tools shipped with this project.
"""
from typing import List, Dict, Any, Set, Optional

_DRUGMATE_CORE_TOOLS = ['terminal',
 'process',
 'read_file',
 'write_file',
 'patch',
 'search_files',
 'vision_analyze',
 'skills_list',
 'skill_view',
 'skill_manage',
 'todo',
 'memory',
 'session_search',
 'clarify',
 'sbdd_confirm',
 'execute_code',
 'delegate_task']

TOOLSETS = {'terminal': {'description': 'Terminal and managed process execution',
              'tools': ['terminal', 'process'],
              'includes': []},
 'file': {'description': 'Read, write, patch and search files',
          'tools': ['read_file', 'write_file', 'patch', 'search_files'],
          'includes': []},
 'vision': {'description': 'Analyze scientific images',
            'tools': ['vision_analyze'],
            'includes': []},
 'skills': {'description': 'Discover and read local skills',
            'tools': ['skills_list', 'skill_view', 'skill_manage'],
            'includes': []},
 'todo': {'description': 'Task planning', 'tools': ['todo'], 'includes': []},
 'memory': {'description': 'Local persistent memory', 'tools': ['memory'], 'includes': []},
 'context_engine': {'description': 'Tools supplied by the active context engine',
                    'tools': [],
                    'includes': []},
 'session_search': {'description': 'Search prior conversations',
                    'tools': ['session_search'],
                    'includes': []},
 'clarify': {'description': 'Ask clarifying questions', 'tools': ['clarify'], 'includes': []},
 'sbdd_confirmation': {'description': 'Review scientific execution decisions',
                       'tools': ['sbdd_confirm'],
                       'includes': []},
 'code_execution': {'description': 'Programmatic tool calls',
                    'tools': ['execute_code'],
                    'includes': []},
 'delegation': {'description': 'Delegate bounded subtasks',
                'tools': ['delegate_task'],
                'includes': []},
 'drugmate-cli': {'description': 'Interactive Python CLI',
                'tools': ['terminal',
                          'process',
                          'read_file',
                          'write_file',
                          'patch',
                          'search_files',
                          'vision_analyze',
                          'skills_list',
                          'skill_view',
                          'skill_manage',
                          'todo',
                          'memory',
                          'session_search',
                          'clarify',
                          'sbdd_confirm',
                          'execute_code',
                          'delegate_task'],
                'includes': []},
 'coding': {'description': 'Workspace tools',
            'tools': ['terminal',
                      'process',
                      'read_file',
                      'write_file',
                      'patch',
                      'search_files',
                      'vision_analyze',
                      'skills_list',
                      'skill_view',
                      'skill_manage',
                      'todo',
                      'memory',
                      'session_search',
                      'clarify',
                      'sbdd_confirm',
                      'execute_code',
                      'delegate_task'],
            'includes': [],
            'posture': True},
 'debugging': {'description': 'Debugging toolkit', 'tools': [], 'includes': ['terminal', 'file']},
 'safe': {'description': 'Image analysis without command execution',
          'tools': [],
          'includes': ['vision']}}

def get_toolset(name: str) -> Optional[Dict[str, Any]]:
    """
    Get a toolset definition by name.
    
    Args:
        name (str): Name of the toolset
        
    Returns:
        Dict: Toolset definition with description, tools, and includes
        None: If toolset not found
    """
    toolset = TOOLSETS.get(name)

    try:
        from tools.registry import registry
    except Exception:
        return toolset if toolset else None

    if toolset:
        merged_tools = sorted(
            set(toolset.get("tools", []))
            | set(registry.get_tool_names_for_toolset(name))
        )
        return {**toolset, "tools": merged_tools}

    return None


def resolve_toolset(name: str, visited: Set[str] = None) -> List[str]:
    """
    Recursively resolve a toolset to get all tool names.
    
    This function handles toolset composition by recursively resolving
    included toolsets and combining all tools.
    
    Args:
        name (str): Name of the toolset to resolve
        visited (Set[str]): Set of already visited toolsets (for cycle detection)
        
    Returns:
        List[str]: List of all tool names in the toolset
    """
    if visited is None:
        visited = set()
    
    # Special aliases that represent all tools across every toolset
    # This ensures future toolsets are automatically included without changes.
    if name in {"all", "*"}:
        all_tools: Set[str] = set()
        for toolset_name in get_toolset_names():
            # Use a fresh visited set per branch to avoid cross-branch contamination
            resolved = resolve_toolset(toolset_name, visited.copy())
            all_tools.update(resolved)
        return sorted(all_tools)

    # Check for cycles / already-resolved (diamond deps).
    # Silently return [] — either this is a diamond (not a bug, tools already
    # collected via another path) or a genuine cycle (safe to skip).
    if name in visited:
        return []

    visited.add(name)

    # Get toolset definition
    toolset = get_toolset(name)
    if not toolset:
        return []

    # Collect direct tools
    tools = set(toolset.get("tools", []))

    # Recursively resolve included toolsets, sharing the visited set across
    # sibling includes so diamond dependencies are only resolved once and
    # cycle warnings don't fire multiple times for the same cycle.
    for included_name in toolset.get("includes", []):
        included_tools = resolve_toolset(included_name, visited)
        tools.update(included_tools)
    
    return sorted(tools)


def resolve_multiple_toolsets(toolset_names: List[str]) -> List[str]:
    """
    Resolve multiple toolsets and combine their tools.
    
    Args:
        toolset_names (List[str]): List of toolset names to resolve
        
    Returns:
        List[str]: Combined list of all tool names (deduplicated)
    """
    all_tools = set()
    
    for name in toolset_names:
        tools = resolve_toolset(name)
        all_tools.update(tools)
    
    return sorted(all_tools)


def get_all_toolsets() -> Dict[str, Dict[str, Any]]:
    """Return definitions for the built-in tool groups."""
    return {name: get_toolset(name) for name in TOOLSETS}


def get_toolset_names() -> List[str]:
    """Return all built-in tool group names."""
    return sorted(TOOLSETS)


def validate_toolset(name: str) -> bool:
    """Return whether a built-in group or all-tools alias is available."""
    return name in {"all", "*"} or name in TOOLSETS


def create_custom_toolset(
    name: str,
    description: str,
    tools: List[str] = None,
    includes: List[str] = None
) -> None:
    """
    Create a custom toolset at runtime.
    
    Args:
        name (str): Name for the new toolset
        description (str): Description of the toolset
        tools (List[str]): Direct tools to include
        includes (List[str]): Other toolsets to include
    """
    TOOLSETS[name] = {
        "description": description,
        "tools": tools or [],
        "includes": includes or []
    }




def get_toolset_info(name: str) -> Dict[str, Any]:
    """
    Get detailed information about a toolset including resolved tools.
    
    Args:
        name (str): Toolset name
        
    Returns:
        Dict: Detailed toolset information
    """
    toolset = get_toolset(name)
    if not toolset:
        return None
    
    resolved_tools = resolve_toolset(name)
    
    return {
        "name": name,
        "description": toolset["description"],
        "direct_tools": toolset["tools"],
        "includes": toolset["includes"],
        "resolved_tools": resolved_tools,
        "tool_count": len(resolved_tools),
        "is_composite": bool(toolset["includes"])
    }




if __name__ == "__main__":
    print("Toolsets System Demo")
    print("=" * 60)
    
    print("\nAvailable Toolsets:")
    print("-" * 40)
    for name, toolset in get_all_toolsets().items():
        info = get_toolset_info(name)
        composite = "[composite]" if info["is_composite"] else "[leaf]"
        print(f"  {composite} {name:20} - {toolset['description']}")
        print(f"     Tools: {len(info['resolved_tools'])} total")
    
    print("\nToolset Resolution Examples:")
    print("-" * 40)
    for name in ["web", "terminal", "safe", "debugging"]:
        tools = resolve_toolset(name)
        print(f"\n  {name}:")
        print(f"    Resolved to {len(tools)} tools: {', '.join(sorted(tools))}")
    
    print("\nMultiple Toolset Resolution:")
    print("-" * 40)
    combined = resolve_multiple_toolsets(["web", "vision", "terminal"])
    print("  Combining ['web', 'vision', 'terminal']:")
    print(f"    Result: {', '.join(sorted(combined))}")
    
    print("\nCustom Toolset Creation:")
    print("-" * 40)
    create_custom_toolset(
        name="my_custom",
        description="My custom toolset for specific tasks",
        tools=["web_search"],
        includes=["terminal", "vision"]
    )
    custom_info = get_toolset_info("my_custom")
    print("  Created 'my_custom' toolset:")
    print(f"    Description: {custom_info['description']}")
    print(f"    Resolved tools: {', '.join(custom_info['resolved_tools'])}")
